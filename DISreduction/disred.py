import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pyraf import iraf
import glob
import os

def test(rawdata_dir='UT160328/',reduction_dir=None,dodetector='blue',dograting='B1200',silent=True):
    print('-'*49)
    print("---------- Drew's DIS reduction script ----------")
    print('-'*49)

    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.twod()
    iraf.apex()
    iraf.longs()
    iraf.oned()

    # check if a reduction directory exists. If not, make one.
    if reduction_dir is None: outdir='reduction_'+rawdata_dir
    check=glob.glob(outdir)
    if len(check)==0:
        cmnd='mkdir '+outdir
        os.system(cmnd)
        print('created reduction directory: '+outdir)

    # establish the read noise and gain values, which are different for blue vs. red
    if dodetector=='blue':
        gain=1.68
        rdnoise=4.9
        label='B'
        badpixfile=''
        fixpix='no'
    if dodetector=='red':
        gain=1.88
        rdnoise=4.6
        label='R'
        badpixfile='badpix_disR.txt'
        fixpix='yes'

    # get an array of the raw image files
    rawfiles = np.array(glob.glob(rawdata_dir+'*fits'))
    nrawfiles = len(rawfiles)

    # get arrays of detector & grating from the headers
    grating = []; detector = []
    for rawfile in rawfiles:
        header = fits.getheader(rawfile)
        grating.append(header['GRATING'])
        detector.append(header['DETECTOR'])
    grating=np.array(grating); detector=np.array(detector)

    # get the desired files based on detector & grating
    gd=np.where((detector==dodetector) & (grating==dograting))
    gdfiles=rawfiles[gd]; ngd=len(gd[0])
    gdfiles.sort()

    # gather some more info from headers
    imtype = []; objname = []; exptime = []; ra = []; dec = []
    for gdfile in gdfiles:
        header = fits.getheader(gdfile)
        imtype.append(header['IMAGETYP'])
        objname.append(header['OBJNAME'])
        exptime.append(header['EXPTIME'])
        ra.append(header['RA'])
        dec.append(header['DEC'])
    imtype=np.array(imtype);   objname=np.array(objname)
    exptime=np.array(exptime); ra=np.array(ra)
    dec=np.array(dec)

    # separate the files into biases, flat, objs, and comps
    bias=np.where(imtype=='zero');  nbias=len(bias[0]); biasfiles=gdfiles[bias]
    flat=np.where(imtype=='flat');  nflat=len(flat[0]); flatfiles=gdfiles[flat]
    obj=np.where(imtype=='object'); nobj=len(obj[0]);   objfiles=gdfiles[obj]
    comp=np.where(imtype=='comp');  ncomp=len(comp[0]); compfiles=gdfiles[comp]

    # Extract the overscan and data sections from the header of whichever file
    biashead=fits.getheader(biasfiles[0])
    biassec=biashead['biassec']
    datasec=biashead['datasec']

    if silent is False:
        # print out some info about the images
        print(str(nbias)+' bias images found')
        for i in range(len(bias[0])): print('   '+biasfiles[i]+',  expt = '+str(exptime[bias][i]))
        print(str(nflat)+' flat images found')
        for i in range(len(flat[0])): print('   '+flatfiles[i]+',  expt = '+str(exptime[flat][i]))

    lamp=[]
    if silent is False: print(str(ncomp)+' comp images found')
    for i in range(len(comp[0])): 
        header = fits.getheader(compfiles[i])
        lamp.append(header['LAMP'])
        if silent is False: print('   '+compfiles[i]+',  expt = '+str(exptime[comp][i])+',  lamp = '+lamp[i])
    lamp=np.array(lamp)

    if silent is False: 
        print(str(nobj)+' object images found')
        for i in range(len(obj[0])): print('   '+objfiles[i]+',  expt = '+str(exptime[obj][i])+
                                           ',  obj = '+objname[obj][i]+',  RA = '+ra[obj][i]+
                                           ',  DEC = '+dec[obj][i])

    ###############################################################################################################
    # average the bias frames using IRAF-IMCOMBINE
    masterzero=outdir+'Zero'+label+'.fits'
    answer=raw_input('step 1. Make master bias (y/n)? ')
    if answer=='y':
        iraf.imcombine(','.join(biasfiles),output=masterzero,combine='average',reject='avsigclip',lsigma='3',
                       hsigma='3',rdnoise=rdnoise,gain=gain)
        print('   master bias made')
    print('-'*50)

    ###############################################################################################################
    # do the overscan and bias subtraction on the flat, comp, and obj frames via IRAF-CCDPROC
    # don't do any trimming for now

    answer=raw_input('step 2. Overscan and bias-subtraction of the flats, comps, and objs. Do it (y/n)? ')
    if answer=='y':
        # set up a list of the input files (flats, comps, objs)
        ff=flatfiles.tolist(); cf=compfiles.tolist(); of=objfiles.tolist()
        infiles=np.array(ff+cf+of)
        infilesfile=outdir+'inlist'+label+'_ccdproc1'
        np.savetxt(infilesfile,infiles,fmt='%s')
        
        # We're making changes to the images now, so need to set up a list of output files
        outfiles=[]
        for i in range(len(infiles)): 
            tmp=infiles[i]; tmp1=tmp.split('/')
            outfiles.append(outdir+'cproc1_'+tmp1[len(tmp1)-1])
        outfiles=np.array(outfiles)
        outfilesfile=outdir+'outlist'+label+'_ccdproc1'
        np.savetxt(outfilesfile,outfiles,fmt='%s')

        # Now for the CCDPROC command
        iraf.ccdproc('@'+infilesfile,output='@'+outfilesfile,ccdtype='',fixpix=fixpix,oversca='yes',trim='yes',
                     zerocor='yes',darkcor='no',flatcor='no',fixfile=badpixfile,biassec=biassec,trimsec=datasec,
                     zero=masterzero,interac='no',low_rej='3',high_re='3')
        print('   overscan and bias subtracting done')
    print('-'*50)

    ###############################################################################################################
    # combine the flats into a master flat
    masterflat=outdir+'Flat'+label+'.fits'
    answer=raw_input('step 3. Make master flat (y/n)? ')
    if answer=='y':
        procfiles=np.array(glob.glob(outdir+'cproc1*fits'))
        imtyp=[]
        for i in range(len(procfiles)):
            head=fits.getheader(procfiles[i])
            imtyp.append(head['imagetyp'])
        imtyp=np.array(imtyp)
        gd=np.where(imtyp=='flat')
        flats=procfiles[gd]

        iraf.imcombine(','.join(flats),output=masterflat,combine='median',reject='avsigclip',lsigma='2',
                       hsigma='2',rdnoise=rdnoise,gain=gain)
    print('-'*50)

    ###############################################################################################################
    # run RESPONSE on the master flat, to make a normalized flat.
    masternormflat=outdir+'NormFlat'+label+'.fits'
    answer=raw_input('step 4. Normalize the response of the master flat (y/n)? ')
    if answer=='y':
        # the normalization spectrum image should probably be set to a subsection of the masternormflat, but i'm using the full image for now.
        iraf.response(calibrat=masterflat,normaliz=masterflat,response=masternormflat,interac='no',functio='spline3',order='5',low_rej='3',high_rej='3')
        print('   normalized master flat created')

    print('-'*50)

    ###############################################################################################################
    # flat field the comps and objs via CCDPROC

    answer=raw_input('step 5. Flat field the comps and objs (y/n)? ')
    if answer=='y':
        # set up a list of the input files (comps, objs)
        tmpfiles=np.array(glob.glob(outdir+'cproc1*fits'))
        imtyp=[]
        for i in range(len(tmpfiles)):
            head=fits.getheader(tmpfiles[i])
            imtyp.append(head['imagetyp'])
        imtyp=np.array(imtyp)
        gd=np.where((imtyp=='object') | (imtyp=='comp'))
        infiles=tmpfiles[gd]
        infilesfile=outdir+'inlist'+label+'_ccdproc2'
        np.savetxt(infilesfile,infiles,fmt='%s')

        # We're making changes to the images now, so need to set up a list of output files
        outfiles=[]
        for i in range(len(infiles)): 
            tmp=infiles[i]; tmp1=tmp.split('/'); tmp2=tmp1[len(tmp1)-1]
            outfiles.append(outdir+tmp2.replace('cproc1','cproc2'))
        outfiles=np.array(outfiles)
        outfilesfile=outdir+'outlist'+label+'_ccdproc2'
        np.savetxt(outfilesfile,outfiles,fmt='%s')

        # Now for the CCDPROC command
        iraf.ccdproc('@'+infilesfile,output='@'+outfilesfile,ccdtype='',fixpix='no',oversca='no',trim='no',
                     zerocor='no',darkcor='no',flatcor='yes',fixfile='',biassec=biassec,trimsec=datasec,
                     flat=masternormflat,interac='no',low_rej='3',high_re='3')
        print('   overscan and bias subtracting done')
    print('-'*50)

    ###############################################################################################################
    # aperture extraction with APALL

    answer=raw_input('step 6. Extract apertures (y/n)? ')
    if answer=='y':
        infiles=np.array(glob.glob(outdir+'cproc2*fits'))
        contpeak=[]
        for infile in infiles:
            head=fits.getheader(infile)
            imdata=fits.getdata(infile)
            print(infile)
            # display the image, and let the user select the continuum peak
            plt.imshow(imdata[1:1028,924:1124], cmap='gray',vmin=-80, vmax=80)
            cont=raw_input('   enter the continuum y pixel value: ')
            print('   the aperture will be centered on line '+str(answer))
#            contpeak.append(answer)



    return answer




