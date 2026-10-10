program main
implicit none
      real::start,finish
      integer,parameter::L=40   !size of the lattice
      real,parameter::s=1.5     !spin number
      real(kind=8)::spin(L,L,L,3),spin2(L,L,L,3)         !spin matrix
      integer::nn1(L,L,L,2,3),nn2(L,L,L,2,3),nn3(L,L,L,2,3)         !nn(nearest neighbor) and nnn
      real,parameter::J1=1.4543D-3,J2=0.6938D-3,J3=-0.11636D-3           !coupling constant, eV
      integer,parameter::MCsteps=20E6,MCequil=5E6,MCsample=1       !MC_steps, MC equillbrium(used in the statistic), and sampling interval
      integer,allocatable::sample(:)                                !MC sampling 
      integer::sampleshape
      real,parameter::Tstart=5,Tend=300,Tsteps=150 !temperature related 
      real(kind=8)::T(Tsteps)
      real(kind=8),parameter::pi=3.141592653589793
      real(kind=8),parameter::kB=8.617333262145D-5 ! boltzmann constant, eV/K
      integer::ii=1,jj=1,kk=1             !iteration number
      integer::a,b,c              !numbers used temporarily
      real(kind=8)::x1,x2,x3,x4
      integer::aa,bb,cc            !lattice index used in MC steps
      real(kind=8)::theta,fhi               !rotate degree of spin
      real(kind=8)::s_temp(3)           !spin temp, used for comparison with the old
      real(kind=8)::dE                  !energy change in each MC step
      real(kind=8)::R                   !acceptance probability
      real(kind=8)::E(0:MCsteps)
      real(kind=8)::mx(MCsteps),my(MCsteps),mz(MCsteps),meanmx(Tsteps),meanmy(Tsteps),meanmz(Tsteps),varmx(Tsteps),varmy(Tsteps),varmz(Tsteps),chix(Tsteps),chiy(Tsteps),chiz(Tsteps),meane(Tsteps),Cv(Tsteps)

call cpu_time(start)
!define nearest neighbor
do ii=1,L,1
    do jj=1,L,1
        do kk=1,L,1

            !first nn1
            if (ii==1) then
                a=ii-1+L
                b=jj
                c=kk
                nn1(ii,jj,kk,1,1)=a
                nn1(ii,jj,kk,1,2)=b
                nn1(ii,jj,kk,1,3)=c
            else
                a=ii-1
                b=jj
                c=kk
                nn1(ii,jj,kk,1,1)=a
                nn1(ii,jj,kk,1,2)=b
                nn1(ii,jj,kk,1,3)=c
            end if

            !second nn1
            if (ii==L) then
                a=ii+1-L
                b=jj
                c=kk
                nn1(ii,jj,kk,2,1)=a
                nn1(ii,jj,kk,2,2)=b
                nn1(ii,jj,kk,2,3)=c
            else
                a=ii+1
                b=jj
                c=kk
                nn1(ii,jj,kk,2,1)=a
                nn1(ii,jj,kk,2,2)=b
                nn1(ii,jj,kk,2,3)=c
            end if

            !first nn2
            if (jj==1) then
                a=ii
                b=jj-1+L
                c=kk
                nn2(ii,jj,kk,1,1)=a
                nn2(ii,jj,kk,1,2)=b
                nn2(ii,jj,kk,1,3)=c
            else
                a=ii
                b=jj-1
                c=kk
                nn2(ii,jj,kk,1,1)=a
                nn2(ii,jj,kk,1,2)=b
                nn2(ii,jj,kk,1,3)=c
            end if


            !second nn2
            if (jj==L) then
                a=ii
                b=jj+1-L
                c=kk
                nn2(ii,jj,kk,2,1)=a
                nn2(ii,jj,kk,2,2)=b
                nn2(ii,jj,kk,2,3)=c
            else
                a=ii
                b=jj+1
                c=kk
                nn2(ii,jj,kk,2,1)=a
                nn2(ii,jj,kk,2,2)=b
                nn2(ii,jj,kk,2,3)=c
            end if

            !first nn3
            if (kk==1) then
                a=ii
                b=jj
                c=kk-1+L
                nn3(ii,jj,kk,1,1)=a
                nn3(ii,jj,kk,1,2)=b
                nn3(ii,jj,kk,1,3)=c
            else
                a=ii
                b=jj
                c=kk-1
                nn3(ii,jj,kk,1,1)=a
                nn3(ii,jj,kk,1,2)=b
                nn3(ii,jj,kk,1,3)=c
            end if
            
            !second nn3
            if (kk==L) then
                a=ii
                b=jj
                c=kk+1-L
                nn3(ii,jj,kk,2,1)=a
                nn3(ii,jj,kk,2,2)=b
                nn3(ii,jj,kk,2,3)=c
            else
                a=ii
                b=jj
                c=kk+1
                nn3(ii,jj,kk,2,1)=a
                nn3(ii,jj,kk,2,2)=b
                nn3(ii,jj,kk,2,3)=c
            end if
        end do
    end do
end do

!initialize temperature
T(1)=Tstart
do ii=2,Tsteps,1
    T(ii)=Tstart+(Tend-Tstart)/(Tsteps-1)*(ii-1)
end do

!initialize spin
spin(:,:,:,1)=0
spin(:,:,:,2)=0
do ii=1,L,1
    if (mod(ii,2)==1) then
        spin(:,:,ii,3)=s
    else
        spin(:,:,ii,3)=-s
    end if
end do
spin2=spin

E(0)=0
!initialize energy
do ii=1,L,1
    do jj=1,L,1
        do kk=1,L,1
             E(0)=E(0)-J1*(spin(ii,jj,kk,1)*(spin(nn1(ii,jj,kk,1,1),nn1(ii,jj,kk,1,2),nn1(ii,jj,kk,1,3),1)+spin(nn1(ii,jj,kk,2,1),nn1(ii,jj,kk,2,2),nn1(ii,jj,kk,2,3),1))+&
                           spin(ii,jj,kk,2)*(spin(nn1(ii,jj,kk,1,1),nn1(ii,jj,kk,1,2),nn1(ii,jj,kk,1,3),2)+spin(nn1(ii,jj,kk,2,1),nn1(ii,jj,kk,2,2),nn1(ii,jj,kk,2,3),2))+&
                           spin(ii,jj,kk,3)*(spin(nn1(ii,jj,kk,1,1),nn1(ii,jj,kk,1,2),nn1(ii,jj,kk,1,3),3)+spin(nn1(ii,jj,kk,2,1),nn1(ii,jj,kk,2,2),nn1(ii,jj,kk,2,3),3)))+&
                      -J2*(spin(ii,jj,kk,1)*(spin(nn2(ii,jj,kk,1,1),nn2(ii,jj,kk,1,2),nn2(ii,jj,kk,1,3),1)+spin(nn2(ii,jj,kk,2,1),nn2(ii,jj,kk,2,2),nn2(ii,jj,kk,2,3),1))+&
                           spin(ii,jj,kk,2)*(spin(nn2(ii,jj,kk,1,1),nn2(ii,jj,kk,1,2),nn2(ii,jj,kk,1,3),2)+spin(nn2(ii,jj,kk,2,1),nn2(ii,jj,kk,2,2),nn2(ii,jj,kk,2,3),2))+&
                           spin(ii,jj,kk,3)*(spin(nn2(ii,jj,kk,1,1),nn2(ii,jj,kk,1,2),nn2(ii,jj,kk,1,3),3)+spin(nn2(ii,jj,kk,2,1),nn2(ii,jj,kk,2,2),nn2(ii,jj,kk,2,3),3)))+&
                      -J3*(spin(ii,jj,kk,1)*(spin(nn3(ii,jj,kk,1,1),nn3(ii,jj,kk,1,2),nn3(ii,jj,kk,1,3),1)+spin(nn3(ii,jj,kk,2,1),nn3(ii,jj,kk,2,2),nn3(ii,jj,kk,2,3),1))+&
                           spin(ii,jj,kk,2)*(spin(nn3(ii,jj,kk,1,1),nn3(ii,jj,kk,1,2),nn3(ii,jj,kk,1,3),2)+spin(nn3(ii,jj,kk,2,1),nn3(ii,jj,kk,2,2),nn3(ii,jj,kk,2,3),2))+&
                           spin(ii,jj,kk,3)*(spin(nn3(ii,jj,kk,1,1),nn3(ii,jj,kk,1,2),nn3(ii,jj,kk,1,3),3)+spin(nn3(ii,jj,kk,2,1),nn3(ii,jj,kk,2,2),nn3(ii,jj,kk,2,3),3)))
        end do
    end do
end do
write(*,*)E(0)

!open file and save data
open(unit=10,file='MT.dat')
write(10,312) '#T','mx','my','mz','chix','chiy','chiz','E','Cv'
312   format(a12' 'a12' 'a12' 'a12' 'a12' 'a12' 'a12' 'a12' 'a12)
311   format(f12.7' 'f12.7' 'f12.7,' ',f12.7,' ',f12.7,' ',f12.7,' ',f12.7,' ',f12.7,' ',f12.7)

!initial sampling serials
sampleshape=MCequil/MCsample+1
allocate(sample(sampleshape))
do ii=1,sampleshape,1
    sample(ii)=MCsteps-MCequil+(ii-1)*MCsample
end do

call random_seed()

!start MC steps
do ii=1,Tsteps,1
    do jj=1,MCsteps,1
        !determine which lattice to rotate
        call random_number(x1)
        call random_number(x2)
        call random_number(x3)
        aa=ceiling(x1*L)
        bb=ceiling(x2*L)
        cc=ceiling(x3*L)
        !determine the rotate angle
        call random_number(theta)
        call random_number(fhi)
        theta=theta*2*pi
        fhi=theta*2*pi
        s_temp(1)=s*sin(theta)*cos(fhi)
        s_temp(2)=s*sin(theta)*sin(fhi)
        s_temp(3)=s*cos(theta)
        !calculate the energy change
        call random_number(x4)
        dE=-J1*(s_temp(1)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),1)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),1))+&
                s_temp(2)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),2)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),2))+&
                s_temp(3)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),3)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),3)))+&
           -J2*(s_temp(1)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),1)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),1))+&
                s_temp(2)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),2)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),2))+&
                s_temp(3)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),3)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),3)))+&
           -J3*(s_temp(1)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),1)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),1))+&
                s_temp(2)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),2)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),2))+&
                s_temp(3)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),3)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),3)))&
           +J1*(spin(aa,bb,cc,1)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),1)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),1))+&
                spin(aa,bb,cc,2)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),2)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),2))+&
                spin(aa,bb,cc,3)*(spin(nn1(aa,bb,cc,1,1),nn1(aa,bb,cc,1,2),nn1(aa,bb,cc,1,3),3)+spin(nn1(aa,bb,cc,2,1),nn1(aa,bb,cc,2,2),nn1(aa,bb,cc,2,3),3)))&
           +J2*(spin(aa,bb,cc,1)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),1)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),1))+&
                spin(aa,bb,cc,2)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),2)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),2))+&
                spin(aa,bb,cc,3)*(spin(nn2(aa,bb,cc,1,1),nn2(aa,bb,cc,1,2),nn2(aa,bb,cc,1,3),3)+spin(nn2(aa,bb,cc,2,1),nn2(aa,bb,cc,2,2),nn2(aa,bb,cc,2,3),3)))&
           +J3*(spin(aa,bb,cc,1)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),1)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),1))+&
                spin(aa,bb,cc,2)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),2)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),2))+&
                spin(aa,bb,cc,3)*(spin(nn3(aa,bb,cc,1,1),nn3(aa,bb,cc,1,2),nn3(aa,bb,cc,1,3),3)+spin(nn3(aa,bb,cc,2,1),nn3(aa,bb,cc,2,2),nn3(aa,bb,cc,2,3),3)))
        R=exp(-dE/(kB*T(ii)))
        if (dE>1 .OR. R > x4) then
            spin(aa,bb,cc,1)=s_temp(1)
            spin(aa,bb,cc,2)=s_temp(2)
            spin(aa,bb,cc,3)=s_temp(3)
            E(jj)=E(jj-1)+dE
        else
            E(jj)=E(jj-1)
        end if
        mx(jj)=sum(spin(:,:,:,1))/(L*L*L)
        my(jj)=sum(spin(:,:,:,2))/(L*L*L) 
        mz(jj)=sum(spin(:,:,:,3))/(L*L*L) 
    end do
    meanmx(ii)=sum(mx(sample))/(sampleshape)
    meanmy(ii)=sum(my(sample))/(sampleshape)
    meanmz(ii)=sum(mz(sample))/(sampleshape)
    varmx(ii)=sum((mx(sample)-meanmx(ii))**2)/(sampleshape)
    varmy(ii)=sum((my(sample)-meanmy(ii))**2)/(sampleshape)
    varmz(ii)=sum((mz(sample)-meanmz(ii))**2)/(sampleshape)
    chix=varmx(ii)/(T(ii))*L*L*L 
    chiy=varmy(ii)/(T(ii))*L*L*L
    chiz=varmz(ii)/(T(ii))*L*L*L
    meane(ii)=sum(E(sample))/(sampleshape)
    Cv(ii)=sum((E(sample)-meane(ii))**2)/(sampleshape*kB*T(ii)*T(ii)*L*L*L)
    write(*,*) 'step=',ii
    write(10,311) T(ii),abs(meanmx(ii)),abs(meanmy(ii)),abs(meanmz(ii)),chix(ii),chiy(ii),chiz(ii),meane(ii),Cv(ii)
    
    !initialize spin
    spin=spin2
end do

!write(*,*) nn(9,2,2,1),nn(9,2,2,2),T(100),T(1),spin(3,3,1),spin(3,3,2),x1,x2
call cpu_time(finish)
write (*,*),'total time',finish-start
end
