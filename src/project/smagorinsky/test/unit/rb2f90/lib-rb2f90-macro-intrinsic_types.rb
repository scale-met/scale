#
#== Fortran90 code generator as macro for intrinsic types of Fortran
#
#Authors::   Yasuhiro MORIKAWA
#Version::   $Id: lib-rb2f90-macro-intrinsic_types.rb,v 1.3 2009-10-12 04:50:16 morikawa Exp $
#Tag Name::  $Name:  $
#Copyright:: Copyright (C) GFD Dennou Club, 2005-2009. All rights reserved.
#License::   See link:../COPYRIGHT
#
#These variables are used to generate f90 code files from Ruby code 
#files. They are expected to help ruby code approximates f90 code for
#as long as possible.
#
#[JAPANESE]
#
#これらの変数は Ruby で記述されたファイルから F90 ファイルを生成
#するための変数です. これらの関数により, できるだけ F90 コードに
#近い形で Ruby のコードを記述できることが期待されます.
#

$type_intent_out = {"Real"    => "real",
                    "Double"  => "real(DP)",
                    "Logical" => "logical",
                    "Int"     => "integer",
                    "Char"    => "character(STRING)"}
$type_intent_in  = {"Real"    => "real",
                    "Double"  => "real(DP)",
                    "Logical" => "logical",
                    "Int"     => "integer",
                    "Char"    => "character(*)"}
$type_intent_inout = $type_intent_in
$type_internal   = {"Real"    => "real",
                    "Double"  => "real(DP)",
                    "Logical" => "logical",
                    "Int"     => "integer",
                    "Char"    => "character(STRING)"}

#
# This format is conformd with "Cprintf" in "dc_string" module
#
$type_fmt        = {"Real"    => "r",
                    "Double"  => "f",
                    "Logical" => "y",
                    "Int"     => "d",
                    "Char"    => "c"}
$type_fmtarg     = {"Real"    => "r",
                    "Double"  => "d",
                    "Logical" => "l",
                    "Int"     => "i",
                    "Char"    => "c1"}

#
# This format is conformd with "SP", "DP" in "dc_types" module
#
$type_numsuf     = {"Real"    => ".0",
                    "Double"  => ".0_DP",
                    "Int"     => ""}

$type_dpsuf     = {"Real"    => "",
                   "Double"  => "_DP"}

#
# This format is conformd with MPI library
#
$mpi_type      = {"Real"    => "MPI_REAL",
                  "Double"  => "MPI_DOUBLE_PRECISION",
                  "Int"     => "MPI_INTEGER", 
                  "Char"    => "MPI_CHARACTER"}
