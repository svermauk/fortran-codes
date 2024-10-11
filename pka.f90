program pka
implicit none
integer :: i,n,ierr
real*8 :: df_model,df_prot,ddf,dpka,pka_model,pka_prot,ratio
character(len=64) :: method,system,prot,protonated,deprotonated,equal,format
real,parameter :: gas_cons = 0.0019872036d0 , temp = 300.d0, ph=7.d0
format = "(A11, A6, 3(F10.1), 3(F8.1), ES12.1, 3X, A18)"

open(unit=1,file="input.dat",status="old")
open(unit=10,file="output.dat",status="unknown")

i=0
do
  read(1,*,iostat=ierr)
  if(ierr/=0)exit
  i=i+1
end do
rewind(1)
n=i-1

write(10,*) "#Method   System   pka(model) dF(model) dF(prot)   ddF    dpka   pka_prot   A-/AH    prot state"

do i=1,n

!  if(i==1) then
!    cycle 
!  end if

  read(1,*) method,system,pka_model,df_model,df_prot
  ddf = df_prot-df_model
  dpka = ddf/(2.303*gas_cons*temp)
  pka_prot = dpka+pka_model

  if(pka_prot.gt.ph) then
    prot="protonated"
  elseif(pka_prot.lt.ph) then
    prot="deprotonated"
  else
    prot=equal
  end if

  ratio = 10.d0**(ph-pka_prot)

  write(10,format) method,system,pka_model,df_model,df_prot,ddf,dpka,pka_prot,ratio,prot
  !write(10,"(11A,6A,F8.1,F8.1,F8.2,F8.2,F8.2,ES15.6,15A)") method,system,pka_model,df_model,df_prot,ddf,dpka,pka_prot,ratio,prot
end do

end program pka
