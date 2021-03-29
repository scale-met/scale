=begin
Congure Converter from version 5.4 to 5.5
=end

def usage
  print "Configure converter from 5.4 to 5.5\n"
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

  # atmos boundary
  if /^&PARAM_ATMOS_BOUNDARY$/i =~ param_name
    print param_name, "\n"
    param_items.each do |item|
      if /ATMOS_BOUNDARY_(START_DATE|UPDATE_DT|INTERP_TYPE)/i !~ item
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
