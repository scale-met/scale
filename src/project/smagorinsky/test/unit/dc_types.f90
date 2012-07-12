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
  ! ���̷��ѥ�᥿���󶡤��ޤ���
  !
  !
  implicit none
  private
  public :: DP
  public :: TOKEN, STRING
  public :: STDIN, STDOUT, STDERR
  public :: NF_KIND_INT1
  public :: NF_KIND_INT2

  integer, parameter:: DP  = kind(0.0d0) ! Double Precision.
                                         ! �����ټ¿����ѿ��μ��̷��ѥ�᥿
                                         ! �Ȥ����Ѥ��ޤ���

  integer, parameter:: TOKEN  = 32       ! Token.
                                         ! ñ��䥭����ɤ��ݻ�����
                                         ! ʸ�����ѿ��μ��̷��ѥ�᥿
                                         ! �Ȥ����Ѥ��ޤ���

  integer, parameter:: STRING = 256      ! String.
                                         ! ʸ������ݻ�����
                                         ! ʸ�����ѿ��μ��̷��ѥ�᥿
                                         ! �Ȥ����Ѥ��ޤ���
                                         !
                                         !--
                                         !��ȯ�Ը�������
                                         !
                                         ! 256 �Ȥ����ͤ˿�����ͳ�Ϥ���ޤ���.
                                         ! ɬ�פʤ�Ф���礭���ͤ�����
                                         ! ���Ƥ⹽���ޤ���.
                                         ! ������ 8 �Х��ȶ����Ȥʤ�褦,
                                         ! 8 ���ܿ��ȤʤäƤ��뤳�Ȥ�
                                         ! �侩���ޤ�.
                                         !
                                         ! SR11000 �κ�Ŭ��
                                         ! FORTRAN90 ����Ѥ������
                                         ! �Ϥ������� 255 �ʲ���
                                         ! ���ꤹ��ɬ�פ�����ޤ�.
                                         !
                                         !++

  integer, parameter:: STDIN  = 5        ! ɸ�����Ϥ������ֹ�
  integer, parameter:: STDOUT = 6        ! ɸ����Ϥ������ֹ�
  integer, parameter:: STDERR = 0        ! ɸ�२�顼���Ϥ������ֹ�

  ! netCDF Fortran ���󥿡��ե������η��ѥ�᥿
  ! (netcdf.inc �ˤ�¸�ߤ��ʤ�)
  !
  integer, parameter:: NF_KIND_INT1 = selected_int_kind(2)
  integer, parameter:: NF_KIND_INT2 = selected_int_kind(4)
end module
