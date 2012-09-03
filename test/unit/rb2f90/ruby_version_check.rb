#!/usr/bin/env ruby
#
#= Ruby version check
#
# Authors::   Yasuhiro MORIKAWA
# Version::   $Id: ruby_version_check.rb,v 1.1 2009-03-25 08:17:36 morikawa Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#Ruby のバージョンチェックを行うためのスクリプトです
#

# このバージョンより大きいか, 等しい場合に終了コード 0 を返す. 
# それ以外の場合は 1 を返す. 
check_ver = "1.8"

script_ver = VERSION

check_ver_ary = check_ver.split('.')
script_ver_ary = script_ver.split('.')

ary_size = check_ver_ary.size
ary_size = script_ver_ary.size if script_ver_ary.size < ary_size

ary_size.times{ |i|
  exit 1 if script_ver_ary[i] < check_ver_ary[i]
}

exit 0
