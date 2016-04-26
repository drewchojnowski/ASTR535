import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pyraf import iraf
import glob
import os
import sys

def test(instr=None,rawdata_dir='UT160327',reduction_dir=None,dograting=None,dofilter=None,silent=False):
    print('='*70)
    print("---------- Drew's APO 3.5m reduction script ----------")
    print('='*70)

    print('loading some basic iraf packages...')
    iraf.noao(); iraf.imred(); iraf.ccdred()
    print('='*70)

##########################################################################################################
    print('1. INSTRUMENT ESTABLISHMENT\n')
    if instr is None:
        print('Which APO 3.5m instrument did you use (default is ARCTIC)? Choices:')
        print(' - enter 1 for ARCTIC (default)')
        print(' - enter 2 for DIS')
        print(' - enter 3 for something else')
        a=raw_input('')
        if a=='':  a='1'
        if a[:1]=='1': instr='arctic'
        if a[:1]=='2': instr='dis'
        if (a[:1]!='1') & (a[:1]!='2'): sys.exit("You're out of luck for now. Sorry!")
        print('\nOk. '+instr.upper()+' data will be reduced')
        print('='*70)

##########################################################################################################
    print('2. RAW DATA DIRECTORY ESTABLISHMENT\n')
    if rawdata_dir is None:
        a=raw_input('Enter the name of the directory where your images are stored: ')
        if a=='': sys.exit('You must specify a raw data directory. Exiting the APO 3.5m reduction script.')
        if a[-1]!='/': a=a+'/'
        rawdata_dir=a
        check=glob.glob(rawdata_dir)
        if len(check)==0: sys.exit('The directory "'+rawdata_dir+'" does not exist. Exiting the APO 3.5m reduction script.')
    else:
        if rawdata_dir[-1:]!='/': rawdata_dir=rawdata_dir+'/'
        check=glob.glob(rawdata_dir)
        if len(check)==0: sys.exit('The directory "'+rawdata_dir+'" does not exist. Exiting the APO 3.5m reduction script.')

    print('Ok. The raw data directory will be "'+rawdata_dir+'"')
    print('='*70)

##########################################################################################################
    print('3. RAW DATA CHECK\n')
    rawfiles = np.array(glob.glob(rawdata_dir+'*fits')); nraw=len(rawfiles)
    if nraw==0: sys.exit('No FITS files found in the raw data directory. Exiting the APO 3.5m reduction script.')
    instrlist=[]; imtype=[]
    for rawfile in rawfiles:
        head=fits.getheader(rawfile)
        instrlist.append(head['instrume'].lower())
    instrlist=np.array(instrlist)
    gd=np.where(instrlist==instr); ngd=len(gd[0])
    if ngd==0: 
        sys.exit('No FITS files from '+instr.upper()+' were found in the raw data directory. Exiting the APO 3.5m reduction script.')
    else:
        rawfiles=rawfiles[gd]
        print('Ok. Found '+str(ngd)+' FITS images with instrument='+instr.upper()+' in the raw data directory.')
        print('='*70)

##########################################################################################################
    print('4. REDUCTION DIRECTORY ESTABLISHMENT\n')
    if reduction_dir is not None:
        check=glob.glob(reduction_dir)
        if len(check)==0: 
            a=raw_input('Specified reduction directory does not exist. Create it (Y/N)? ')
            if (a=='') | (a[:1].lower()=='y'): 
                cmnd='mkdir '+reduction_dir
                os.system(cmnd)
                print('"'+reduction_dir+'" has been created.')
    else:
        reduction_dir=rawdata_dir[:-1]+'_reduction/'
        print('The reduction directory will be '+reduction_dir)
        print(' - hit ENTER to accept')
        print(' - otherwise enter a different reduction directory name: ')
        a=raw_input('')
        if a=='':
            check=glob.glob(reduction_dir)
            if len(check)==0:
                cmnd='mkdir '+reduction_dir
                os.system(cmnd)
                print('\nOk. "'+reduction_dir+'" has been created.')
        else:
            reduction_dir=a
            check=glob.glob(reduction_dir)
            if len(check)==0:
                cmnd='mkdir '+reduction_dir 
                os.system(cmnd)
                print('\nOk. "'+reduction_dir+'" has been created.')

    print('\nOk. Reduction products will be stored in '+reduction_dir+'.')
    print('='*70)

##########################################################################################################
    if instr=='dis':
        print('loading some more IRAF packages for DIS reduction...')
        iraf.twod(); iraf.apex(); iraf.longs(); iraf.oned()
        print('='*70)
        print('5. DIS PARAMETER ESTABLISHMENT\n')
        grating=[]; imagetypes=[]
        for rawfile in rawfiles:
            head=fits.getheader(rawfile)
            grating.append(head['grating'])
            imagetypes.append(head['imagetyp'].lower())
        imagetypes=np.array(imagetypes)
        grating=np.array(grating)
        uniquegrating=np.unique(grating)

        if dograting is None:
            print("Which DIS grating do you want to reduce? Here's what you have:")
            for i in range(len(uniquegrating)): 
                if i==0: 
                    print(' - enter '+str(i+1)+' for '+uniquegrating[i]+' (default)')
                else:
                    print(' - enter '+str(i+1)+' for '+uniquegrating[i])
            a=raw_input('')
            if a=='': a='1'
            dograting=uniquegrating[int(a[:1])-1]
            gd=np.where(grating==dograting); ngd=len(gd[0])
            rawfiles=rawfiles[gd]; imagetypes=imagetypes[gd]
            print('\nOk. '+dograting+' data will be reduced: '+str(ngd)+' total images')
        else:
            gd=np.where(grating==dograting); ngd=len(gd[0])
            if ngd==0:
                sys.exit('No '+dograting+' images found. Exiting the APO 3.5m reduction script.')
            else:
                rawfiles=rawfiles[gd]; imagetypes=imagetypes[gd]
                print('\nOk. '+dograting+' images will be reduced: '+str(ngd)+' total images')

        bias=np.where((imagetypes=='bias') | (imagetypes=='zero')); nbias=len(bias[0])
        zerocor='yes'
        if nbias==0:
            print('\nNo bias images were found. Proceed without bias correction? (y/n)')
            a=raw_input('')
            if (a[:1]=='y') | (a==''): 
                zerocor='no'; biasfiles=''; nbias=0
            else:
                sys.exit('Exiting the APO 3.5m reduction script.')
        else:
            biasfiles=rawfiles[bias]; nbias=len(bias[0])

        # get values depending on R vs B detector
        detector=dograting[:1].lower()
        if detector=='b':
            gain=1.68
            rdnoise=4.9
            badpixfile=''
            fixpix='no'
        if detector=='r':
            gain=1.88
            rdnoise=4.6
            badpixfile='badpix_disR.txt'
            fixpix='yes'

        head=fits.getheader(rawfiles[0])
        biassec=head['biassec']
        datasec=head['datasec']

        genericlabel=dograting

        print('\nOk. DIS parameters established.')
        print('='*70)

##########################################################################################################
    if instr=='arctic':
        print('5. ARCTIC PARAMETER ESTABLISHMENT\n')
        gain=2.0
        rdnoise=3.7
        badpixfile=''
        fixpix='no'

        filters=[]; imagetypes=[]
        for rawfile in rawfiles:
            head=fits.getheader(rawfile)
            filters.append(head['filter'])
            imagetypes.append(head['imagetyp'].lower())
        imagetypes=np.array(imagetypes)
        filters=np.array(filters)
        uniquefilters=np.unique(filters)

        bias=np.where(imagetypes=='bias'); nbias=len(bias[0])
        zcor='yes'
        if nbias==0:
            print('\nNo bias images were found. Proceed without bias correction? (y/n)')
            a=raw_input('')
            if (a[:1]=='y') | (a==''): 
                zcor='no'; biasfiles=''
            else:
                sys.exit('Exiting the APO 3.5m reduction script.')
        else:
            biasfiles=rawfiles[bias]; nbias=len(bias[0])

        # establish ARCTIC filter
        if dofilter is None:
            print('Which ARCTIC filter do you want to reduce? ')
            for i in range(len(uniquefilters)): 
                if i==0: 
                    print(' - enter '+str(i+1)+' for '+uniquefilters[i]+' (default)')
                else:
                    print(' - enter '+str(i+1)+' for '+uniquefilters[i])
            a=raw_input('')
            if a=='': a='1'
            dofilter=uniquefilters[int(a[:1])-1]
            gd=np.where(filters==dofilter); ngd=len(gd[0])
            rawfiles=rawfiles[gd]
            print('\nOk. '+dofilter+' filter images will be reduced: '+str(ngd)+' total images')
        else:
            gd=np.where(filters==dofilter); ngd=len(gd[0])
            if ngd==0:
                sys.exit('No '+dofilter+' images found. Exiting the APO 3.5m reduction script.')
            else:
                rawfiles=rawfiles[gd]
                print('\nOk. '+dofilter+' filter images will be reduced: '+str(ngd)+' total images')


        head=fits.getheader(rawfiles[0])
        biassec=head['bsec11']
        datasec=head['dsec11']

        genericlabel=dofilter.replace(' ','_')
        genericlabel=genericlabel.replace('/','_')

        print('\nOk. ARCTIC parameters established.')
        print('='*70)

##########################################################################################################
    print('6. GET IMAGE TYPES\n')
    nrawfiles = len(rawfiles)

    # gather some more info from headers
    imtype = []; exptime = []; ra = []; dec = []
    for rawfile in rawfiles:
        header = fits.getheader(rawfile)
        imtype.append(header['IMAGETYP'].lower())
        exptime.append(header['EXPTIME'])
        ra.append(header['RA'])
        dec.append(header['DEC'])
    imtype=np.array(imtype)
    exptime=np.array(exptime); ra=np.array(ra)
    dec=np.array(dec)

    # separate the files into biases, flat, objs, and comps
    flat=np.where(imtype=='flat');  nflat=len(flat[0]); flatfiles=rawfiles[flat]
    obj=np.where(imtype=='object'); nobj=len(obj[0]);   objfiles=rawfiles[obj]
    if instr=='dis':
        comp=np.where(imtype=='comp'); ncomp=len(comp[0]); compfiles=rawfiles[comp]
        print('Ok. You have '+str(nbias)+' biases, '+str(nflat)+' flats, '+str(ncomp)+' comps, and '+str(nobj)+' objects.')
    if instr=='arctic':
        print('Ok. You have '+str(nbias)+' biases, '+str(nflat)+' flats, and '+str(nobj)+' objects.')
    print('='*70)


##########################################################################################################
    print('7. MAKE MASTER BIAS\n')
    if zcor=='yes':
        masterzero=reduction_dir+'Zero_'+genericlabel+'.fits'
        a=raw_input('Average combine biases into master bias (y/n)? ')
        if (a[:1]=='y') | (a==''):
            iraf.imcombine(','.join(biasfiles),output=masterzero,combine='average',reject='avsigclip',lsigma='3',
                           hsigma='3',rdnoise=rdnoise,gain=gain)
            print('\nOk. A master bias has been made: '+masterzero)
        print('='*70)
    else:
        masterzero=''
        print('\nOk, you have no bias files. Proceeding')


##########################################################################################################
    print('8. BIAS AND OVERSCAN CORRECTION\n')
    a=raw_input('Overscan and bias-subtraction of the flats, objs, and if DIS, comps. Do it (y/n)? ')
    if (a[:1]=='y') | (a==''):
        # set up a list of the input files (flats, comps, objs)
        if instr=='dis': infiles=np.concatenate((flatfiles,compfiles,objfiles),axis=0)
        if instr=='arctic': infiles=np.concatenate((flatfiles,objfiles),axis=0)

        infilesfile=reduction_dir+'inlist_'+genericlabel+'_ccdproc1'
        np.savetxt(infilesfile,infiles,fmt='%s')
        
        # We're making changes to the images now, so need to set up a list of output files
        outfiles=[]
        for infile in infiles:
            tmp=infile.split('/')
            tmp1=tmp[len(tmp)-1].split('.fits')
            outfiles.append(reduction_dir+tmp1[0]+'_cproc1.fits')
        outfiles=np.array(outfiles)
        outfilesfile=reduction_dir+'outlist_'+genericlabel+'_ccdproc1'
        np.savetxt(outfilesfile,outfiles,fmt='%s')

        iraf.ccdproc('@'+infilesfile,output='@'+outfilesfile,ccdtype='',fixpix=fixpix,oversca='yes',trim='yes',
                     zerocor=zcor,darkcor='no',flatcor='no',fixfile=badpixfile,biassec=biassec,trimsec=datasec,
                     zero=masterzero,interac='no',low_rej='3',high_re='3')
        print('\nOk. Overscan and bias subtracting done.')
    print('='*70)

##########################################################################################################
    print('9. MAKE MASTER FLAT\n')
    masterflat=reduction_dir+'Flat_'+genericlabel+'.fits'
    a=raw_input('Median combine the flats (y/n)? ')
    if (a[:1]=='y') | (a==''):
        procfiles=np.array(glob.glob(reduction_dir+'*_cproc1.fits'))
        imtyp=[]
        for i in range(len(procfiles)):
            head=fits.getheader(procfiles[i])
            imtyp.append(head['imagetyp'].lower())
        imtyp=np.array(imtyp)
        gd=np.where(imtyp=='flat')
        flats=procfiles[gd]

        iraf.imcombine(','.join(flats),output=masterflat,combine='median',reject='avsigclip',lsigma='2',
                       hsigma='2',rdnoise=rdnoise,gain=gain)
        print('\nOk. A master flat has been made: '+masterflat)
    print('='*70)

##########################################################################################################
    print('10. MAKE NORMALIZED MASTER FLAT\n')
    masternormflat=reduction_dir+'NormFlat_'+genericlabel+'.fits'
    a=raw_input('Normalize the response of the master flat (y/n)? ')
    if (a[:1]=='y') | (a==''):
        if instr=='arctic':
            im=fits.getdata(masterflat)
            imdim=len(im)
            meanval=np.mean(im[imdim-100:imdim+100,imdim-100:imdim+100])
            iraf.imarith(operand1=masterflat,op='/',operand2=meanval,result=masternormflat)
            print('\nOk. A normalized master flat has been created: '+masternormflat)
        if instr=='dis':
            iraf.response(calibrat=masterflat,normaliz=masterflat+'[*,400:500]',response=masternormflat,
                          interac='no',functio='spline3',order='5',low_rej='3',high_rej='3')
            print('\nOk. A normalized master flat has been created: '+masternormflat)
    print('='*70)

##########################################################################################################
    print('11. FLAT FIELD THE IMAGES\n')

    a=raw_input('Flat field the objs & comps by the normalized master flat (y/n)? ')
    if (a[:1]=='y') | (a==''):
        tmpfiles=np.array(glob.glob(reduction_dir+'*_cproc1.fits'))
        imtyp=[]
        for i in range(len(tmpfiles)):
            head=fits.getheader(tmpfiles[i])
            imtyp.append(head['imagetyp'].lower())
        imtyp=np.array(imtyp)
        gd1=np.where(imtyp=='object')
        if instr=='dis':
            gd2=np.where(imtyp=='comp')
            infiles=np.concatenate((tmpfiles[gd1],tmpfiles[gd2]),axis=0)
        if instr=='arctic':
            infiles=tmpfiles[gd1]
        infilesfile=reduction_dir+'inlist_'+genericlabel+'_ccdproc2'
        np.savetxt(infilesfile,infiles,fmt='%s')

        outfiles=[]
        for infile in infiles: 
            outfiles.append(infile.replace('cproc1','cproc2'))
        outfiles=np.array(outfiles)
        outfilesfile=reduction_dir+'outlist_'+genericlabel+'_ccdproc2'
        np.savetxt(outfilesfile,outfiles,fmt='%s')

        iraf.ccdproc('@'+infilesfile,output='@'+outfilesfile,ccdtype='',fixpix='no',oversca='no',trim='no',
                     zerocor='no',darkcor='no',flatcor='yes',fixfile='',biassec=biassec,trimsec=datasec,
                     flat=masternormflat,interac='no',low_rej='3',high_re='3')
        print('\nOk. The objs and comps have been flat fielded.')
    print('='*70)


##########################################################################################################
    if instr=='dis':
        print('Next steps (done manually for now):')
        print('='*70)
        print('12. Extract apertures via APALL.')
        print('13. Run IDENTIFY on one or more comps.')
        print('14. (optional) Run REIDENTIFY to transfer solution among comps.')
        print('15. Assign solution(s) to objs using REFSPECTRA.')
        print('16. Run DISPCOR on the objs to apply the wavelength solution.')
        print('='*70)
        print('Other things that should be done:')
        print('='*70)
        print('step ?1. flux calibration.')
        print('step ?2. sky/background subtraction.')
        print('step ?3. heliocentric velocity correction.')


    return imtyp



