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

  if /^\s*&PARAM_(.+)$/ =~ line
    params.push [line,[]]
    inflag = $1.strip.upcase
  elsif /^\s*&.+\/\s*$/ =~ line
    params.push [line]
  elsif /^\s*\// =~ line
    inflag = false
  elsif inflag
    params[-1][1].push line
  else
    params.push line
  end

  if /ATMOS_PHY_TB_TYPE\s*=\s*["'](\w+)['"]/i =~ line
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
    if /ATMOS_PHY_TB_HYBRID_SGS_TYPE\s*=\s*["'](.+)['"]/ =~ line
      tb = $1
    else
      tb = "SMAGORINSKY"
    end
    if /ATMOS_PHY_TB_HYBRID_PBL_TYPE\s*=\s*["'](.+)['"]/ =~ line
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

  param_name  = param[0]
  param_items = param[1]

  # HISTORY
  ## scale_file_history (old: gtool_history)
  if /&PARAM_HISTORY/ =~ param_name
    print "&PARAM_FILE_HISTORY\n"
    param_items.each do |item|
      print item.sub(/^(\s*)HISTORY_/, '\1FILE_HISTORY_'), "\n"
    end
    print "/\n"
    next
  end
  if /&HISTITEM(.*\s+)item=(.+)$/ =~ param_name
    print "&HISTORY_ITEM#{$1}name=#{$2}\n"
    next
  end
  if /&PARAM_HIST\s*$/ =~ param_name
    print "&PARAM_FILE_HISTORY_CARTESC\n"
    param_items.each do |item|
      item.sub!("HIST_BND", "FILE_HISTORY_CARTESC_BOUNDARY")
      print item.sub("HIST_", "FILE_HISTORY_CARTESC"), "\n"
    end
    print "/\n"
    next
  end

  # GRID
  ## scale_atmos_grid_cartesC_index (old: scale_grid_index)
  param_name.sub!("&PARAM_INDEX", "&PARAM_ATMOS_GRID_CARTESC_INDEX")
  ## scale_atmos_grid_cartesC (old: scale_grid_cartesian)
  param_name.sub!("&PARAM_GRID", "&PARAM_ATMOS_GRID_CARTESC")
  ["OCEAN", "LAND", "URBAN"].each do |model|
    ## scale_(model)_grid_cartesC_index (old: scale_(model)_grid_index)
    param_name.sub!("&PARAM_#{model}_INDEX", "&PARAM_#{model}_GRID_CARTESC_INDEX")
    ## scale_(model)_grid_cartesC (old: scale_(model)_grid_cartesian)
    param_name.sub!("&PARAM_#{model}_GRID", "&PARAM_#{model}_GRID_CARTESC")
  end

  # MAP Projection
  ## scale_mapprojection (old: scale_mapproj)
  if /&PARAM_MAPPROJ/ =~ param_name
    print "&PARAM_MAPPROJECTION\n"
    param_items.each do |item|
      print item.sub(/^(\s*)MPRJ_/, '\1MAPPROJECTION_'), "\n"
    end
    print "/\n"
    next
  end

  # Metrics
  if /&PARAM_GTRANS/ =~ param_name
    print "&PARAM_ATMOS_GRID_CARTESC_METRIC\n"
    param_items.each do |item|
      print item.sub(/^(\s*)GTRANS_/, '\1ATMOS_GRID_CARTESC_METRIC_'), "\n"
    end
    print "/\n"
    next
  end

  # Dynamics
  ## mod_atmos_dyn
  if /&PARAM_ATMOS_DYN/ =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /^(\s*)ATMOS_DYN_enable_coriolis\s*=\s*([^,\s]+)/i =~ item && $2 == ".true."
        print "#{$1}ATMOS_DYN_coriolis_type = \"SPHERE\",\n"
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end

  # Boundary-Layer scheme
  if /&PARAM_ATMOS\s*$/ =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /^(\s*)ATMOS_PHY_TB_TYPE\s*=\s*/ =~ item
        print "#{$1}ATMOS_PHY_TB_TYPE = \"#{tb}\",\n" if tb
        print "#{$1}ATMOS_PHY_BL_TYPE = \"#{bl}\",\n" if bl
      else
        print item, "\n"
      end
    end
    print "/\n"
    next
  end
  if /&PARAM_TIME/ =~ param_name
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
  if /&PARAM_ATMOS_PHY_TB_HYBRID/ =~ param_name
    next
  end

  # Nesting
  if /&PARAM_NEST/ =~ param_name
    print "&PARAM_COMM_CARTESC_NEST\n"
    param_items.each do |item|
      if /^(\s*)USE_NESTING\s*=/ !~ item && /^(\s*)OFFLINE\s*=/ !~ item
        print item, "\n"
      end
    end
    print "/\n"
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
