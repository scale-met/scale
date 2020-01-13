=begin
Congure Converter from version 5.3 to 5.4
=end

def usage
  print "Configure converter from 5.3 to 5.4\n"
  print "Usage: ruby #$0 ord.conf\n"
  print "Reuslt will be output to the standard output\n"
  exit
end

fname = ARGV.shift


params = Array.new
inflag = false
File.foreach(fname) do |line|
  line.chomp!

  if /^\s*&PARAM_(.+)$/i =~ line # namelist name
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

end



# print parameters after conversion
params.each do |param|

  if String === param
    print param, "\n"
    next
  end

  param_name  = param[0].strip
  param_items = param[1]

  # coriolis parameter
  if /^&PARAM_ATMOS_DYN$/i =~ param_name
    coriolis = Array.new
    print "&PARAM_ATMOS_DYN\n"
    param_items.each do |item|
      if /coriolis/i =~ item
        coriolis.push item.sub("ATMOS_DYN_","")
      else
        print item, "\n"
      end
    end
    print "/\n"
    if coriolis.any?
      print "\n&PARAM_CORIOLIS\n"
      coriolis.each do |item|
        print item, "\n"
      end
      print "/\n"
    end
    next
  end


  # topo
  if /^&PARAM_TOPO$/i =~ param_name
    print "&PARAM_TOPOGRAPHY\n"
    param_items.each do |item|
      print item.sub("TOPO", "TOPOGRAPHY"), "\n"
    end
    print "/\n"
    next
  end


  # KUSAKA01
  if /^&PARAM_URBAN_DYN_KUSAKA01$/i =~ param_name
    print param_name, "\n"
    list = Array.new
    param_items.each do |item|
      if /DTS_MAX/i =~ item || /BOUND/i =~ item
        print item, "\n"
      else
        list.push item
      end
    end
    if list.any?
      dir = File.dirname(fname)
      fn = File.join(dir, "param.kusaka01.dat")
      if File.exist?(fn)
        $stderr.print "[ERROR] 'param.kusaka01.dat' file exists\n"
        exit(-1)
      end
      File.open(fn, "w") do |file|
        file.print "&PARAM_URBAN_DATA\n"
        list.each do |item|
          item.sub!(/AH(L?)/, 'AH\1_TBL')
          file.print item, "\n"
        end
        file.print "/\n"
      end
      print " URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME = '#{File.basename(fn)}',\n"
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
