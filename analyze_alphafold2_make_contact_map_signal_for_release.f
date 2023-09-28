
C simple code that measures all atomic contacts between a signal and the
C rest of the structure

C the code has five arguments:

C (1) alphafold2 .pdb name - contains full-length structure
C (2) distance threshold for atomic contacts (we used 4A)
C (3) number of residues in the signal
C (4) pLDDT score threshold for counting contacts
C (5) number of residues either side of cut-site to omit from calcs

C output is written to:

C (1) signal_contacts.txt

      program analyze_alphafold

      implicit real (a-h,o-z)
      common/aherandom/idumahe
      integer, allocatable :: icontact(:,:)
      integer, allocatable :: jcontact(:,:)
      integer, allocatable :: kcontact(:,:)
      integer ninres(10000)
      real bfac(10000)
      real x1(10000,100) ! max 10000 residues, 100 atoms in each
      real y1(10000,100)
      real z1(10000,100)
      character*3 jnk,knk,nam        
      character*1 lnk,lnk_lst
      character*50 junk
      character*50 file1,file1pseudo
      character*30 initial_pdb
      character*50 file2,file2pseudo
      character*50 filefinal
      character*80 string   
      parameter                                                     
     &         (pi=3.141592654,                                        
     &       twopi=6.283185307,                                       
     &         eps=0.001,                                               
     &        eps2=0.999999,                                           
     &        eps3=0.0001,                                             
     &      dtorad=0.0174532925)   
      call getarg(1,junk)
      read(junk,*) file1     
      call getarg(2,junk)
      read(junk,*) cut_dist 
      call getarg(3,junk)
      read(junk,*) nsignal
      call getarg(4,junk)
      read(junk,*) bfac_thresh ! only count contacts when both residues
C                                have bfacs>= this value
      call getarg(5,junk)
      read(junk,*) nskip ! skip these res before and after cutsite
     

      cut_dist2=cut_dist**2

      nframe=1 ! AHE changed here for Venkat to use AHEMODEl structure
      natm=0
      do j=1,10000 
        ninres(j)=0
      enddo
          
C read in the pdb file to identify all the residues

      do n=1,30
        initial_pdb(n:n)=' '
      enddo

      nres_sig_acc=0
      nres_rem_acc=0

      len1=len(trim(file1)) ! WTF this is overkill
ccc   write(*,*)'working on ',file_pdb(1:len1)
ccc   write(*,*)'working on ',initial_pdb(1:len1)

      initial_pdb(1:len1)=file1(1:len1)

      open(unit=10,file=file1,status='unknown')
10    read(10,12,end=15)string
12    format(a80)
      if(string(1:4).eq.'ATOM') then    
        read(string,14)ires,x,y,z,occ,bfac_tmp
        natm=natm+1
        ninres(ires)=ninres(ires)+1
14      format(22x,i4,4x,3f8.3,2f6.2)
        bfac(ires)=bfac_tmp
        if(ninres(ires).eq.1) then
          if(ires.le.nsignal) then
            nres_sig_acc=nres_sig_acc+1
          else
            nres_rem_acc=nres_rem_acc+1
          endif
        endif
        x1(ires,ninres(ires))=x
        y1(ires,ninres(ires))=y
        z1(ires,ninres(ires))=z
        goto 10
      else
        goto 10
      endif
15    close(10)

c     do i=1,ires
c       write(*,*)'check ',i,bfac(i)
c     enddo

      write(*,*)' just read ',natm,' atoms  and ',ires,
     &              ' residues'
       
! do the analysis

      ncont_atm_sig_rest=0
      ncont_res_sig_rest=0

C loop over all residues in the signal up to the cutsite-nskip

      do i1=1,nsignal-nskip

C skip this residue if its pLDDT score isn't good enough

        if(bfac(i1).lt.bfac_thresh) cycle

C loop over all residues in the protein after the cutsite+nskip

        do j1=nsignal+nskip+1,ires ! ignore i:i+1,i:i+2

          ididthisres=0

C skip this residue if its pLDDT score isn't good enough

          if(bfac(j1).lt.bfac_thresh) cycle

C loop over all atoms in the current signal residue

          do i2=1,ninres(i1)

C loop over all atoms in the current body residue

            do j2=1,ninres(j1)

C measure the (squared) distance

              dist2=(x1(i1,i2)-x1(j1,j2))**2+
     &              (y1(i1,i2)-y1(j1,j2))**2+
     &              (z1(i1,i2)-z1(j1,j2))**2

              if(dist2.le.cut_dist2) then
                ncont_atm_sig_rest=ncont_atm_sig_rest+1
                ididthisres=1
              endif

            enddo
          enddo

          if(ididthisres.eq.1) ncont_res_sig_rest=
     &                         ncont_res_sig_rest+1

        enddo

      enddo

      open(unit=11,file='signal_contacts.txt',status='unknown')

      write(11,801)initial_pdb(1:30),nsignal,ires,
     &             nres_sig_acc,nres_rem_acc,
     &             ncont_res_sig_rest,ncont_atm_sig_rest
801   format('pdbname ',a30,' #res_sig: ',i6,' #res_tot: ',i6,
     &       ' #res_sig_acc: ',i6,' #res_rem_acc: ',i6,
     &       ' #res-res-conts: ',i6,' #atm-atm-conts: ',i6)

      close(11)

      stop
      end
