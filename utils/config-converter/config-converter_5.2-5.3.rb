=begin
Congure Converter from version 5.2 to 5.3
=end

def usage
  print "Configure converter from 5.2 to 5.3\n"
  print "Usage: ruby #$0 ord.conf\n"
  print "Reuslt will be output to the standard output\n"
  exit
end

fname = ARGV.shift


# check turbulent scheme
hybrid = false
tb = nil
bl = nil
params = Array.new
inflag = false
File.foreach(fname) do |line|
  line.chomp!

  if /^\s*&PARAM_(.+)$/i =~ line # namelist name
    params.push [line,[]]
    inflag = $1.strip.upcase
  elsif /^\s*&NM_MP_SN14_(.+)$/i =~ line # namelist name (irregular)
    params.push [line,[]]
    inflag = $1.strip.upcase
  elsif /^\s*&.+\/\s*$/ =~ line # namelist oneline
    params.push [line]
  elsif /^\s*\// =~ line # namelist endline
    inflag = false
  elsif inflag # namelist item
    params[-1][1].push line
  else # other line
    params.push line
  end

  if /ATMOS_PHY_TB_TYPE\s*=\s*[""](\w+)[""]/i =~ line
    case $1.upcase
    when "HYBRID"
      hybrid = true
    when "MYNN"
      bl = $1
      hybrid = false
    else
      tb = $1
      hybrid = false
    end
  end

  if hybrid
    if /ATMOS_PHY_TB_HYBRID_SGS_TYPE\s*=\s*[""](.+)[""]/i =~ line
      tb = $1
    else
      tb = "SMAGORINSKY"
    end
    if /ATMOS_PHY_TB_HYBRID_PBL_TYPE\s*=\s*[""](.+)[""]/i =~ line
      bl = $1
    else
      bl = "MYNN"
    end
  end

end



# print parameters after conversion
params.each do |param|

  if String === param
    print param, "\n"
    next
  end

  param_name  = param[0].strip
  param_items = param[1]

  # PRC & COMM
  if /^&PARAM_PRC$/i =~ param_name
    print "&PARAM_PRC_CARTESC\n"
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_COMM$/i =~ param_name
    print "&PARAM_COMM_CARTESC\n"
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_NEST$/i =~ param_name
    print "&PARAM_COMM_CARTESC_NEST\n"
    param_items.each do |item|
      if /^(\s*)USE_NESTING\s*=/i !~ item && /^(\s*)OFFLINE\s*=/i !~ item
        print item, "\n"
      end
    end
    print "/\n"
    next
  end

  # FILE, HISTORY, MONIT, EXTIN
  if /^&PARAM_FILEIO$/i =~ param_name
    print "&PARAM_FILE_CARTESC\n"
    param_items.each do |item|
      print item.sub(/FILEIO/i, "FILE_CARTESC"), "\n"
    end
    print "/\n"
    next
  end
  ## scale_file_history (old: gtool_history)
  if /^&PARAM_HISTORY$/i =~ param_name
    print "&PARAM_FILE_HISTORY\n"
    param_items.each do |item|
      print item.sub(/HISTORY_/i, "FILE_HISTORY_"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_HIST$/i =~ param_name
    print "&PARAM_FILE_HISTORY_CARTESC\n"
    param_items.each do |item|
      item.sub!(/HIST_BND/i, "FILE_HISTORY_CARTESC_BOUNDARY")
      print item.sub(/HIST_/i, "FILE_HISTORY_CARTESC_"), "\n"
    end
    print "/\n"
    next
  end
  if /&HISTITEM(.*\s+)ITEM=(.+)$/i =~ param_name
    print "&HISTORY_ITEM#{$1}name=#{$2}\n"
    next
  end
  if /&MONITITEM(.*\s+)ITEM=(.+)$/i =~ param_name
    print "&MONITOR_ITEM#{$1}name=#{$2}\n"
    next
  end
  if /&EXTITEM(.+)$/i =~ param_name
    print "&EXTERNAL_ITEM#{$1}\n"
    next
  end

  # GRID
  ## scale_atmos_grid_cartesC_index (old: scale_grid_index)
  if /^&PARAM_INDEX$/i =~ param_name
    param_name = "&PARAM_ATMOS_GRID_CARTESC_INDEX"
  end
  ## scale_atmos_grid_cartesC (old: scale_grid_cartesian)
  if /^&PARAM_GRID$/i =~ param_name
    param_name = "&PARAM_ATMOS_GRID_CARTESC"
  end
  ["OCEAN", "LAND", "URBAN"].each do |model|
    ## scale_(model)_grid_cartesC_index (old: scale_(model)_grid_index)
    if /^&PARAM_#{model}_INDEX$/i =~ param_name
      param_name = "&PARAM_#{model}_GRID_CARTESC_INDEX"
    end
    ## scale_(model)_grid_cartesC (old: scale_(model)_grid_cartesian)
    if /^&PARAM_#{model}_GRID$/i =~ param_name
      param_name = "&PARAM_#{model}_GRID_CARTESC"
    end
  end

  # Tracer
  if /^&PARAM_TRACER$/i =~ param_name
    next
  end
  if /^&PARAM_TRACER_KAJINO13$/i =~ param_name
    next
  end

  # MAP Projection
  ## scale_mapprojection (old: scale_mapproj)
  if /^&PARAM_MAPPROJ$/i =~ param_name
    print "&PARAM_MAPPROJECTION\n"
    param_items.each do |item|
      print item.sub(/^(\s*)MPRJ_/i, "\1MAPPROJECTION_"), "\n"
    end
    print "/\n"
    next
  end

  # Metrics
  if /^&PARAM_GTRANS$/i =~ param_name
    print "&PARAM_ATMOS_GRID_CARTESC_METRIC\n"
    param_items.each do |item|
      print item.sub(/^(\s*)GTRANS_/i, "\1ATMOS_GRID_CARTESC_METRIC_"), "\n"
    end
    print "/\n"
    next
  end

  # Dynamics
  ## mod_atmos_dyn
  if /^&PARAM_ATMOS_DYN$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /^(\s*)ATMOS_DYN_ENABLE_CORIOLIS\s*=\s*([^,\s]+)/i =~ item && /\.true\./i =~ $2
        print "#{$1}ATMOS_DYN_CORIOLIS_TYPE = \"SPHERE\",\n"
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end

  # Cloud Microphysics scheme
  if /^&PARAM_BIN$/i =~ param_name
    print "&PARAM_ATMOS_PHY_MP_SUZUKI10_bin\n"
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
    next
  end
  if /^&NM_MP_SN14(.*)$/i =~ param_name
    print param_name.sub(/NM_MP_SN14_/i, "PARAM_ATMOS_PHY_MP_SN14_"), "\n"
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_ATMOS_PHY_MP_BIN2BULK$/i =~ param_name
    next
  end

  # Surface flux scheme

  if /^&PARAM_ATMOS_PHY_SF$/i =~ param_name
    print "&PARAM_ATMOS_PHY_SF_BULK\n"
    param_items.each do |item|
      print item.sub(/ATMOS_PHY_SF_/i, "ATMOS_PHY_SF_BULK_"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_ATMOS_PHY_SF_BULKCOEF$/i =~ param_name
    next
  end

  # Boundary-Layer scheme
  if /^&PARAM_ATMOS$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /^(\s*)ATMOS_PHY_TB_TYPE\s*=\s*/i =~ item
        print "#{$1}ATMOS_PHY_TB_TYPE = \"#{tb}\",\n" if tb
        print "#{$1}ATMOS_PHY_BL_TYPE = \"#{bl}\",\n" if bl
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end
  if /^&PARAM_TIME$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /^(\s*)TIME_DT_ATMOS_PHY_TB(\s*=.+)$/ =~ item
        print "#{$1}TIME_DT_ATMOS_PHY_TB#{$2}\n" if tb
        print "#{$1}TIME_DT_ATMOS_PHY_BL#{$2}\n" if bl
      elsif /^(\s*)TIME_DT_ATMOS_PHY_TB_UNIT(\s*=.+)$/ =~ item
        print "#{$1}TIME_DT_ATMOS_PHY_TB_UNIT#{$2}\n" if tb
        print "#{$1}TIME_DT_ATMOS_PHY_BL_UNIT#{$2}\n" if bl
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end
  if /^&PARAM_ATMOS_PHY_TB_HYBRID$/i =~ param_name
    next
  end
  if /^&PARAM_ATMOS_PHY_TB_MYNN$/i =~ param_name
    next
  end

  # Land
  if /^&PARAM_LAND$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      next if /LAND_DO/i =~ item
      if /LAND_TYPE\s*=\s*[""](THIN-)?SLAB[""]/i =~ item
        print " LAND_DYN_TYPE = \"BUCKET\",\n"
        print " LAND_SFC_TYPE = \"SKIN\",\n"
      elsif /LAND_TYPE\s*=\s*[""]THICK-SLAB[""]/i =~ item
        print " LAND_DYN_TYPE = \"BUCKET\",\n"
        print " LAND_SFC_TYPE = \"COPY\",\n"
      elsif /LAND_TYPE\s*=\s*[""]CONST[""]/i =~ item
        print " LAND_DYN_TYPE = \"CONST\",\n"
        print " LAND_SFC_TYPE = \"COPY\",\n"
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end
  if /^&PARAM_LAND_PHY_SLAB$/i =~ param_name
    print "&PARAM_LAND_DYN_BUCKET\n"
    param_items.each do |item|
      print item.sub(/PHY_SLAB/i, "DYN_BUCKET").sub(/PHY_UPDATE/i, "DYN_BUCKET_UPDATE"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_LAND_SFC_SLAB$/i =~ param_name
    print "&PARAM_LAND_SFC_SKIN\n"
    param_items.each do |item|
      print item.sub(/SLAB/i, "SKIN"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_LAND_SFC_THIN_SLAB$/i =~ param_name
    print "&PARAM_LAND_SFC_SKIN\n"
    param_items.each do |item|
      print item.sub(/THIN_SLAB/i, "SKIN"), "\n"
    end
    print "/\n"
    next
  end
  if /&PARAM_LAND_SFC_THICK_SLAB$/i =~ param_name
    print "&PARAM_LAND_SFC_FIXED_TEMP\n"
    param_items.each do |item|
      print item.sub(/THICK_SLAB/i, "FIXED_TEMP"), "\n"
    end
    print "/\n"
    next
  end

  # Ocean
  if /^&PARAM_OCEAN$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      next if /OCEAN_DO/i =~ item
      print item.sub(/OCEAN_TYPE/i, "OCEAN_DYN_TYPE"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_OCEAN_PHY_SLAB$/i =~ param_name
    print "&PARAM_OCEAN_DYN_SLAB\n"
    param_items.each do |item|
      print item.sub(/PHY/i, "DYN"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_OCEAN_PHY_FILE$/i =~ param_name
    print "&PARAM_OCEAN_DYN_SLAB\n"
    param_items.each do |item|
      print item.sub(/PHY/i, "DYN").sub(/FILE/i, "SLAB"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_ROUGHNESS$/i =~ param_name
    print "&PARAM_OCEAN_ROUGHNESS\n"
    param_items.each do |item|
      print item.sub(/ROUGHNESS_TYPE/i, "OCEAN_RGN_TYPE"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_ROUGHNESS_MILLER92$/i =~ param_name
    print "&PARAM_OCEAN_PHY_ROUGHNESS_MILLER92\n"
    param_items.each do |item|
      print item.sub(/ROUGHNESS/i, "OCEAN_PHY_ROUGHNESS"), "\n"
    end
    print "/\n"
    next
  end
  if /^&PARAM_ROUGHNESS_MOON07$/i =~ param_name
    print "&PARAM_OCEAN_PHY_ROUGHNESS_MOON07\n"
    param_items.each do |item|
      print item.sub(/ROUGHNESS/i, "OCEAN_PHY_ROUGHNESS"), "\n"
    end
    print "/\n"
    next
  end

  # Urban
  if /^&PARAM_URBAN$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      next if /URBAN_DO/i =~ item
      if /URBAN_TYPE\s*=\s*[""]SLC[""]/i =~ item
        print " URBAN_DYN_TYPE = \"KUSAKA01\",\n"
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end
  if /^&PARAM_URBAN_PHY_SLC$/i =~ param_name
    print "&PARAM_URBAN_DYN_KUSAKA01\n"
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
    next
  end

  # Restart
  if /^&PARAM_RESTART$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      next if /RESTART_RUN/i =~ item
      print item, "\n"
    end
    print "/\n"
    next
  end

  # Mkinit
  if /^&PARAM_SBMAERO$/i =~ param_name
    next
  end
  if /^&PARAM_MKINIT_INTERPORATION$/i =~ param_name
    next
  end

  # Others
  print param_name, "\n"
  if param_items
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
  end

end
