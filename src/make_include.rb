include Math


if ARGV.length < 4
  print <<EOL
Usage: ruby #$0 DXYZ LX LY LZ [LZB] [IBLOCK] [JBLOCK]
       DXYZ  : resolution (m)
       LX    : domain length of x-direction (m)
       LY    : domain length of y-direction (m)
       LZ    : domain length of z-direction (m)
       LZB   : depth of buffer region (m)
       IBLOCK: block size of x-direction
       JBLOCK: block size of y-direction
EOL
  exit 0
end

dxyz = ARGV.shift.to_f
lx   = ARGV.shift.to_f
ly   = ARGV.shift.to_f
lz   = ARGV.shift.to_f
lzb  = ( ARGV.shift || lz/2 ).to_f
iblock  = ARGV.shift
jblock  = ARGV.shift

imax = (lx/dxyz).ceil
jmax = (ly/dxyz).ceil

k1 = (lz/dxyz).ceil
fact = 1.1
k2 = (log(1.0 - lzb*(1-fact)/dxyz)/log(fact)).ceil
kmax = k1 + k2

iblock ||= imax
jblock ||= jmax

nhalo = 2

class Float
  alias :__to_s :to_s
  def to_s
    s = self.__to_s
    /e/ =~ s ? s.sub(/e/,"E") : s+"D0"
  end
end

print <<EOS
  !-----------------------------------------------------------------------------
  !
!++ scale3 grid parameters (#{dxyz}m res., #{lz+lzb}km model top)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO = #{nhalo} ! # of halo cells: z
  integer, private, parameter :: IHALO = #{nhalo} ! # of halo cells: x
  integer, private, parameter :: JHALO = #{nhalo} ! # of halo cells: y

  real(8), private, parameter :: DX  = #{dxyz} ! length in the main region [m]: x
  real(8), private, parameter :: DY  = #{dxyz} ! length in the main region [m]: y
  real(8), private, parameter :: DZ  = #{dxyz} ! length in the main region [m]: z

  real(8), private, parameter :: BUFFER_DZ =  #{lzb} ! thickness of buffer region [m]: z
  real(8), private, parameter :: BUFFER_DX =  0.0D0 ! thickness of buffer region [m]: x
  real(8), private, parameter :: BUFFER_DY =  0.0D0 ! thickness of buffer region [m]: y
  real(8), private, parameter :: BUFFFACT  =  #{fact} ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = #{kmax} ! # of computational cells: z
  integer, private, parameter :: IMAX = #{imax} ! # of computational cells: x
  integer, private, parameter :: JMAX = #{jmax} ! # of computational cells: y

  integer, private, parameter :: KA   = #{kmax+2*nhalo} ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   = #{imax+2*nhalo} ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   = #{jmax+2*nhalo} ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   = #{nhalo+1} ! start point of inner domain: z, local
  integer, private, parameter :: KE   = #{kmax+nhalo} ! end   point of inner domain: z, local
  integer, private, parameter :: IS   = #{nhalo+1} ! start point of inner domain: x, local
  integer, private, parameter :: IE   = #{imax+nhalo} ! end   point of inner domain: x, local
  integer, private, parameter :: JS   = #{nhalo+1} ! start point of inner domain: y, local
  integer, private, parameter :: JE   = #{jmax+nhalo} ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  = #{imax*jmax} ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  =  1 !
  integer, private, parameter :: IJE  = #{imax*jmax} !

  integer, private, parameter :: IBLOCK = #{iblock} !
  integer, private, parameter :: JBLOCK = #{jblock} !
EOS
