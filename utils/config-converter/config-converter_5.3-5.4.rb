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


  # Others
  print param_name, "\n"
  if param_items
    param_items.each do |item|
      print item, "\n"
    end
    print "/\n"
  end

end
