
      implicit none
      integer,parameter :: ni=1000
      real(8) :: ax,asi(3,ni),aa(3,3),bb(3,3),Va,Pi,aai(3,3)
      real(8) :: al(3),al1(3),rsi(3,ni)
      integer :: i,MKI,MI,nmkd,mkd(ni),Kion(ni)

      Pi=4.d0*atan(1.d0)

      read(*,*) ax
      read(*,*) aa(:,1)
      read(*,*) aa(:,2)
      read(*,*) aa(:,3)
      read(*,*) MKI,MI,nmkd
      do i=1,MI
         read(*,*) Kion(i),rsi(:,i),mkd(i)
      end do

      aa=ax*aa

      al(1)=sqrt(sum(aa(:,1)**2))
      al(2)=sqrt(sum(aa(:,2)**2))
      al(3)=sqrt(sum(aa(:,3)**2))
      write(*,*) "a1",al(1)
      write(*,*) "a2",al(2)
      write(*,*) "a3",al(3)

      al1(1)=10.d0
      al1(2)=10.d0
      al1(3)=10.d0
      do i=1,3
         aa(:,i)=aa(:,i)*al1(i)/al(i)
      end do
      al(1)=sqrt(sum(aa(:,1)**2))
      al(2)=sqrt(sum(aa(:,2)**2))
      al(3)=sqrt(sum(aa(:,3)**2))
      write(*,*) "a1",al(1)
      write(*,*) "a2",al(2)
      write(*,*) "a3",al(3)

      Va=aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1)
     &  +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1)
     &  -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)

      bb(1,1) = aa(2,2)*aa(3,3) - aa(3,2)*aa(2,3) ! x-component
      bb(2,1) = aa(3,2)*aa(1,3) - aa(1,2)*aa(3,3) ! y
      bb(3,1) = aa(1,2)*aa(2,3) - aa(2,2)*aa(1,3) ! z
      bb(1,2) = aa(2,3)*aa(3,1) - aa(3,3)*aa(2,1)
      bb(2,2) = aa(3,3)*aa(1,1) - aa(1,3)*aa(3,1)
      bb(3,2) = aa(1,3)*aa(2,1) - aa(2,3)*aa(1,1)
      bb(1,3) = aa(2,1)*aa(3,2) - aa(3,1)*aa(2,2)
      bb(2,3) = aa(3,1)*aa(1,2) - aa(1,1)*aa(3,2)
      bb(3,3) = aa(1,1)*aa(2,2) - aa(2,1)*aa(1,2)
      bb(:,:)=bb(:,:)*2.d0*Pi/Va

      aai(:,:)=transpose(bb(:,:))/(2.d0*Pi)

      asi(1:3,1:MI)=matmul( aai(1:3,1:3),rsi(1:3,1:MI) )

      ax=maxval(al1(1:3))
      aa=aa/ax

      write(*,*) ax
      write(*,'(1x,3f20.14)') aa(:,1)
      write(*,'(1x,3f20.14)') aa(:,2)
      write(*,'(1x,3f20.14)') aa(:,3)
      write(*,*) MKI,MI,nmkd
      do i=1,MI
         write(*,'(1x,i3,2x,3f18.12,i3)') Kion(i),asi(:,i),mkd(i)
      end do

      aa=ax*aa
      rsi(:,1:MI)=matmul(aa,asi(:,1:MI) )

      do i=1,MI
         write(*,*) rsi(:,i)
      end do

      stop
      end


