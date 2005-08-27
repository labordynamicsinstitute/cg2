! grouping program
! See description at end of file
! requires module my_data (my_data.lars.f90)
      program groups
      use MY_DATA
      implicit none
      CHARACTER(LEN=*), PARAMETER :: this_stamp = "Time-stamp: <01/03/19 creec003>"
      CHARACTER(LEN=*), PARAMETER :: orig_author = "Rob Creecy"
      CHARACTER(LEN=*), PARAMETER :: this_author = "Lars Vilhuber, Rob Creecy"
      CHARACTER(LEN=*), PARAMETER :: this_version = "0.3"

      INTEGER*4 i,j, thispers, thisfirm, person, afirm , nextfirm, err

      write(6,01) orig_author, this_author, this_version, this_stamp
01    format("Grouping program, original author:",A,/,&
             "Modified by: ",A,/,"Version: ",A,/,A)

      read(5,*) version,ncells,npers,nfirm,ncov
      write(6,10) ncells,npers,nfirm
10    format("Number of distinct cells ",I10,/, &
             "Number of persons      ",I10,/, &
             "Number of firms        ",I10,/)
      
      ALLOCATE(byp(2,ncells),byf(2,ncells),m(npers+nfirm),mtype(npers+nfirm), &
                pg(npers), fg(nfirm), pindex(npers), findex(nfirm), stat = err)
      if (err .eq. 0) write(6,20)
20    format(" Allocated integer variables ")
      ALLOCATE(ptraced(npers),ftraced(nfirm) , ponstack(npers), fonstack(nfirm),stat=err)
      if (err .eq. 0) write(6,30)
30    format(" Allocated logical variables ")
      ptraced = .FALSE.
      ftraced = .FALSE.
      ponstack = .FALSE.
      fonstack = .FALSE.
      m = 0
      mtype = 0
      write(6,40)
40    format(" Initialized variables")
      read(7,*) (byp(1,i),byp(2,i), i = 1,ncells)
      read(8,*) (byf(1,i),byf(2,i), i = 1,ncells)
      write(6,50)
 50   format(" Finished reading data, starting preprocessing")
      do j = 1,ncells
        pindex(byp(1,j)) = j
        findex(byf(2,j)) = j
      end do
      g = 1
      call check_data()
! Start with firm 1
      nextfirm = 1
      mpoint = 1
      m(mpoint) = 1
      mtype(mpoint) = 2
      fonstack(1) = .TRUE.
      do while (mpoint .gt. 0)
!         if (g .gt. 1) write(6,*) g, mpoint, m(mpoint),mtype(mpoint)
         call trace()
         if (mpoint .eq. 0) then
            g = g + 1
            do while (nextfirm .lt. nfirm .and. fg(nextfirm) .ne. 0)
               nextfirm = nextfirm + 1
            end do
            if (fg(nextfirm) .eq. 0) then
               mpoint = 1
               m(mpoint) = nextfirm
               mtype(mpoint) = 2
               fonstack(nextfirm) = .TRUE.
            end if
         end if
      end do
      write(6,60)
60    format(" Finished preprocessing, writing data.")
      do j = 1,ncells
         write(4,*) byp(1,j),byp(2,j),pg(byp(1,j))
      end do
      write(6,70)
70    format(" Finished writing groups to Unit 4")
      stop

      end

      subroutine check_data()
	use MY_DATA
      implicit none
      INTEGER*4 i,j, thispers, thisfirm, person, afirm, ncellcount
!     check persons are sorted 1-npers, and ncells is correct
      ncellcount = 0
      thispers = 1
      thisfirm = 1
      do i = 1,ncells
        if (byp(1,i) .ne. thispers) then
	  if (byp(1,i) .ne. thispers + 1) then
	    write(6, 80) thispers, byp(1,i),i
80	    format("Error - by person file not sorted or missing sequence number",/, &
                   "Previous person # - ",i10,"  This person # - ",i10,"  index in file - ",i10)
          stop
          end if
          thispers = thispers + 1
        end if
        if (byf(2,i) .ne. thisfirm) then
	  if (byf(2,i) .ne. thisfirm + 1) then
	    write(6, 90) thisfirm, byf(2,i),i
90	    format("Error - by firm file not sorted or missing sequence number",/, &
                   "Previous firm # - ",i10,"  This firm # - ",i10,"  index in file - ",i10)
          stop
          end if
          thisfirm = thisfirm + 1
        end if
      end do
      write(6,100)
100   format("Data Checked - By Person and By Firm files correctly sorted and sequenced")
      return
      end

      subroutine trace()
        use MY_DATA
      implicit none
      INTEGER*4 i,j, thispers, thisfirm, person, afirm, lower, upper
        if (mtype(mpoint) .eq. 2) then
           thisfirm = m(mpoint)
           mpoint = mpoint - 1
           fg(thisfirm) = g
           ftraced(thisfirm) = .TRUE.
           fonstack(thisfirm) = .FALSE.
           if (thisfirm .eq. 1) then
             lower = 1
           else
              lower = findex(thisfirm - 1) + 1
           end if
           upper = findex(thisfirm)
           do person = lower,upper
              thispers = byf(1,person)
              pg(thispers) = g
              if (.not. ptraced(thispers) .and. .not. ponstack(thispers)) then
                 mpoint = mpoint + 1
                 m(mpoint) = thispers
                 mtype(mpoint) = 1
                 ponstack(thispers) = .TRUE.        
              end if
           end do

        else
           thispers = m(mpoint)
           mpoint = mpoint - 1
           pg(thispers) = g
           ptraced(thispers) = .TRUE.
           ponstack(thispers) = .FALSE.
           if (thispers .eq. 1) then
              lower = 1
           else
              lower = pindex(thispers - 1) + 1
           end if
           upper = pindex(thispers)                      
           do afirm = lower,upper
              thisfirm = byp(2,afirm)
              fg(thisfirm) = g
              if (.not. ftraced(thisfirm) .and. .not. fonstack(thisfirm)) then
                 mpoint = mpoint + 1
                 m(mpoint) = thisfirm
                 mtype(mpoint) = 2
                 fonstack(thisfirm) = .TRUE.
              end if
           end do

        end if
        return
        end
   
! ===== BEGIN DESCRIPTION ===== 
! Taken from groupdoc.txt by Rob Creecy
! Description of the groups program.
! 
! An issue in linear modeling of person and firm effects is how
! to determine which effects in the model are estimable. 
! 
! The groups program helps get at this goal by dividing a
! file of observations on persons and firms into distinct groups.
! A group is defined as all persons and firms that are
! connected through some migration of persons between
! firms in that group, and such that there is no migration
! of a person within the group to any firm outside the group.
! 
! A recursive definition for group 1 is:
!   
!   Firm 1 is in group g = 1
!   Repeat until no more persons or firms are added
!     Add all persons employed by a firm in group g to group g
!     Add all firms that have employed a person in group g to group g
!     
! Other groups are found in the same way, starting with the first
! firm not already assigned to a group.
! 
! The claim is that all the person and firm effects within a group
! except one, which can be arbitrarily set to zero, are estimable.
! (A proof would be nice)
! 
! Detailed description of the code.
! 
! The essential input is a list of distinct person-firm cell indicies
! that occur within the data file to be modeled. e.g.
! 
!        Person    Firm
!           1       2678
!           1       4163
!           1      10226
!           2        961
!           2       2511
!           2       2881
!           2       6856
!           2      10226
!           .        .
!           .        .
! 
! This data is stored in a file sorted by person (e.g. cellsbypers.txt)
! and in a separate file sorted by firm (e.g. cellsbyfirm.txt)
! 
! These files are read into the arrays byp(2,ncells) and byf(2,ncells)
! respectively, so that the first row of each array contains person
! indicies and the second row contains firm indicies.
! 
! The problem size is specified by:
! 
!   ncells npers nfirm
! 
! These three numbers are read on stdin (FORTRAN unit 5)
! and used to allocate arrays of the correct size. They
! may be stored in a file, groups.in, and redirected
! to stdin when the program is invoked.
!  
!   ncells is the number of distinct person/firm cells, 
!   and should be equal to the number of records in the
!   cellsbyxxx input files.
! 
!   npers is the number of persons, indicies from 1 to npers
!   nfirm is the number of firms, indicies from 1 to nfirm
! 
! Output is a file just like the cells by person input file,
! but with a third column containing the assigned group number.
! It is written on FORTRAN unit 4, and is conventionally called
! groupsbyp.txt .
! 
! The program is run as follows (in c shell):
! 
! % setenv FORT7 cellsbyp.txt
! % setenv FORT8 cellsbyf.txt
! % setenv FORT4 groupsbyp.txt      
! % ./groups < groups.in > groups.log
!  
! Implemented algorithm
! 
! The recursive definition of groups above needs to be expanded
! with more detail in order to have a working, efficient algorithm.
! A method needs to be developed to keep track of which persons and firms
! have been assigned to a group and also which persons and firms remain
! to be examined ("traced") to see if they have additional firms or
! persons which need to be added to the group. A person group array (pg)
! and a firm group array (fg) keep track of the groups assigned to persons and
! firms. A group of zero means the group has yet to be assigned.
! 
! The definition above requires that we be able to easily find "all persons
! employed by a firm" and "all firms that have employed a person." The
! by person (byp) and by firm (byf) arrays contain in sequence just the
! needed information. These arrays need to be indexed with a person index
! (pindex indexes byp) or a firm index (findex indexes byf) that can be used
! to find the first and last persons or firms in a sequence. These indices
! are computed so they point to the last person or firm in the sequence.
! 
! A stack is used to keep track of which persons and firms have been identified
! as part of the current group but still need to be "traced" to see if any
! additional firms or persons need to be added to the group. The stack has
! two components, the array m which stores a firm or person index that
! needs to be traced, and the array mtype which indicates if the corresponding
! index is a person index (mtype = 1) or firm index (mtype = 2).  A stack pointer,
! mpoint, indexes the top of the stack. When mpoint is zero the stack is empty
! and the group is complete.
! 
! Two additional logical arrays ponstack and fonstack help improve the program's efficiency
! by indicating if a person index or firm index is currently on the stack. The
! logical arrays ptraced and ftraced indicate if a person or firm index has
! already been traced, helping to prevent repeated tracing of the same person
! or firm.
! 
! groups() {
!    group g = 1
!    push firm 1 onto the stack
!    while the stack is not empty {
!      remove the person or firm from the top of the stack
!      trace() that person or firm, adding firms or persons to the group and stack
!      if the stack is empty {
!        group g = g + 1
!        find the next firm not assigned to a group and push it on the stack
!        if all firms have been assigned, finished
!       }
!    }
! }
!        
! trace() {
!      if this is a firm to trace {
!         assign this firm to group g
!         mark this firm  as traced
!         find all persons in this firm (using findex)
!         for each person found {
!           assign this person to group g
!           if this person has not been traced and it is not already on the stack,
!             push this person onto the stack
!         }
!      }
!      else this is a person to trace { 
!         assign this person to group g
!         mark this person as traced
!         find all firms for this person (using pindex)
!         for each firm found {
!           assign this firm to group g
!           if this firm has not been traced and it is not already on the stack,
!             push this firm onto the stack
!         }
!      }
! }
! 
