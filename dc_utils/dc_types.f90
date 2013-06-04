!== Kind type parameter value
!
! Authors::   Yasuhiro MORIKAWA, Eizi TOYODA
! Version::   $Id: dc_types.f90,v 1.1 2009-03-20 09:09:52 morikawa Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2000-2005. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
! This file provides dc_types
!

module dc_types
  !
  !== Overview
  !
  ! 種別型パラメタを提供します。
  !
  !
  implicit none
  private
  public :: DP, SP
  public :: TOKEN, STRING
  public :: STDIN, STDOUT, STDERR
  public :: NF_KIND_INT1
  public :: NF_KIND_INT2

  integer, parameter:: DP  = kind(0.0d0) ! Double Precision.
                                         ! 倍精度実数型変数の種別型パラメタ
                                         ! として用います。
  integer, parameter:: SP  = kind(0.0e0) ! Single Precision.

  integer, parameter:: TOKEN  = 32       ! Token.
                                         ! 単語やキーワードを保持する
                                         ! 文字型変数の種別型パラメタ
                                         ! として用います。

  integer, parameter:: STRING = 256      ! String.
                                         ! 文字列を保持する
                                         ! 文字型変数の種別型パラメタ
                                         ! として用います。
                                         !
                                         !--
                                         !開発者向け情報
                                         !
                                         ! 256 という値に深い理由はありません.
                                         ! 必要ならばより大きな値を設定
                                         ! しても構いません.
                                         ! ただし 8 バイト境界となるよう,
                                         ! 8 の倍数となっていることを
                                         ! 推奨します.
                                         !
                                         ! SR11000 の最適化
                                         ! FORTRAN90 を使用する場合に
                                         ! はだいたい 255 以下に
                                         ! 指定する必要があります.
                                         !
                                         !++

  integer, parameter:: STDIN  = 5        ! 標準入力の装置番号
  integer, parameter:: STDOUT = 6        ! 標準出力の装置番号
  integer, parameter:: STDERR = 0        ! 標準エラー出力の装置番号

  ! netCDF Fortran インターフェイスの型パラメタ
  ! (netcdf.inc には存在しない)
  !
  integer, parameter:: NF_KIND_INT1 = selected_int_kind(2)
  integer, parameter:: NF_KIND_INT2 = selected_int_kind(4)
end module
