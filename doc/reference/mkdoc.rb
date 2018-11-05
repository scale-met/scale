=begin
Document generator for SCALE

This script generate list of namelist group and history.
The document is output in the HTML format.

Usage: ruby mkdoc.rb doc_dir
=end
require "pp"

Encoding.default_external = Encoding::UTF_8

def usage
  print "Usage: #$0 doc_dir\n"
  exit
end

output_dir = "./dox"

srcdir = ARGV.shift || usage
title = ARGV.shift || File.basename(File.expand_path(File.join(`pwd`,"..")))


def parse_line(line)
  line_org = line.dup
  line.sub!(/^\s*&/, "") # ignore '&' at top of line

  # get comment
  if /["'][^"']*![^"']*["']/ !~ line && /^([^!]*)!(.*)$/ =~ line
    line = $1
    comment = $2
    comment.sub!(/^\s*</,"") if comment
    comment.strip! if comment
  else
    comment = nil
  end

  if /^(.+)&\s*$/ =~ line # '&' at end of line indicates continuation line
    cont = true
    reg = $1
  else
    cont = false
    reg = line
    if /&/ =~ line && /["']([^"']*&[^"']*)+["']/ !~ line
      # if '&' is found (except in strings), it indicates a parsing error.
      p line_org
      p line
      raise "parse error"
    end
  end
  reg.strip!
  return [comment, cont, reg]
end

def parse_dir(dir)
  files = Array.new
  Dir[File.join(dir,"*")].each do |d|
    next if /\A\.{1,2}\Z/ =~ d
    if File.directory?(d)
      files << parse_dir(d)
    elsif /\.(F|f)90\Z/ =~ d
      files.push d
    end
  end
  files
end

def get_data(name, data, vars)
  if /\A([^\(]+)\(([^\)]+)\)\Z/ =~ name
    var = $1
    idx = $2
    if ary = data[var]
      if idx.to_i.to_s == idx
        return ary[id.to_i] || name
      else
        if (h = vars[idx.upcase]) && h[:type]=="integer"
          return ary[h[:val].to_i-1] || name
        end
      end
    end
  end
  return name
end

Index = {
  "I_DENS" => {:type=>"integer", :val=>1},
  "I_MOMZ" => {:type=>"integer", :val=>2},
  "I_MOMX" => {:type=>"integer", :val=>3},
  "I_MOMY" => {:type=>"integer", :val=>4},
  "I_RHOT" => {:type=>"integer", :val=>5},
  "I_QTRC" => {:type=>"integer", :val=>6},
}

files = parse_dir(srcdir)

files = files.flatten.sort!

tree = Hash.new
nm_params = Hash.new
history = Hash.new
files.each do |fname|

  mod_f = fname.sub(/#{srcdir}\//,"").sub(/\.[Ff]90\Z/,"")

  vars = Hash.new
  namelists = Hash.new
  hist = Hash.new
  modname = nil
  fnames = [fname]
  data = Hash.new
  if /scale_(ae)?tracer_(.+)\.F90\Z/ =~ fname && $2 != "sdm"
    fnames.push File.join(File.dirname(fname),"../../../include","inc_#{$1}tracer_#{$2}.h")
  end
  fnames.each do |fn|
    File.open(fn) do |file|
      while line = file.gets
        cont = true
        comments = Array.new
        line_new = ""
        while cont
          com, cont, li = parse_line(line)
          line_new << li
          comments.push(com) if com
          line = file.gets if cont
        end
        case line_new
        when /^\s*(public|private)\s* ::/i
          next
        when /^\s*(?:program|module)\s*([^\s]+)\s*$/i
          modname = $1.strip
        when /^([^,]+).*(intent\([^\)]+\))?.*::([^=]+)(=.*)?$/i
          next if $2
          type = $1
          name = $3
          val = $4
          type.strip!
          name = name.strip.upcase
          val.sub!(/\A=/,"").strip! if val
          while /\A([^,(]+)(?:\(([^)]+)\))?,?(.*)\Z/ =~ name
            na = $1.strip
            dim = $2 && $2.strip
            name = $3 && $3.strip
            ty = type
            ty = ty + ", dimension(#{dim})" if dim
            vars[na] = {:type => ty, :val => nil, :comments => comments}
          end
          vars[na][:val] = val if val
        when /NAMELIST\s*\/([^\)]+)\/(.+)$/i
          group = $1
          lists = $2
          group.strip!
          lists = lists.split(",").map{|c| c.strip.upcase}
          namelists[group] = lists
        when /call FILE_HISTORY_in\s*\((.+)$/i
          str = $1.strip.sub(/\)\Z/,"").strip
          next if /scale_file_history/ =~ modname
          str.gsub!(/\(:[^)]*\)/,'')
          info = str.split(",").map{|c| c.strip.sub(/\A'?(.*)'\Z/,'\1')}
          hist[info[1]] = {:unit => info[3], :desc => info[2], :var => info[0]}
        when /call FILE_HISTORY_reg\s*\((.+)$/i
          str = $1.strip.sub(/\)\Z/,"").strip
          next if /scale_file_history/ =~ modname
          str.gsub!(/\(:[^)]*\)/,'')
          info = str.split(",").map{|c| c.strip.sub(/\A'?(.*)'\Z/,'\1')}
          hist[info[0]] = {:unit => info[2], :desc => info[1], :var => info[0]}
        when /data\s+([^\s]+)\s+\/(.+)\/$/i
          name = $1
          body = $2
          if /Vinfo/ =~ body
            ary = Array.new
            while body
              if /\A(Vinfo\(.+?\)),\s*?(Vinfo.*)\Z/ =~ body
                ary.push $1
                body = $2.strip
              elsif /\A(Vinfo\(.+\))\Z/ =~ body
                ary.push $1
                body = nil
              else
                raise
              end
            end
            data[name] = ary
          else
            data[name] = body.split(',').map{|s| s.strip.sub(/^'/,"").sub(/'$/,"")}
          end
        end

      end # while line
    end # open
  end

  vars.update(Index){|k,v1,v2| v1}

  if namelists.empty? && hist.empty?
    # warn "#{File.basename(fname)} namelist or history not found: namelist #{namelists.any?}, history #{hist.any?}"
    next
  end

  parent = File.dirname(mod_f)
  system("mkdir -p #{output_dir}/#{parent}")

  tree[parent] ||= Array.new
  tree[parent].push File.basename(mod_f)

  File.open("#{output_dir}/#{mod_f}.dox","w") do |file|

    file.print <<EOL
!> @par NAMELIST
!>    <ul>
EOL
    if namelists.empty?
      file.print "!>      <li>No namelist group</li>\n"
    end
    namelists.each do |group, lists|
      file.print <<EOL
!>      <li id="#{group}">#{group}</li>
!> @anchor namelist_#{modname}_#{group}
!>      <table border=1>
!>        <tr><th>name</th><th>type</th><th>default value</th><th>comment</th></tr>
EOL
      lists.each do |v|
        unless info = vars[v]
          info = nil
          vars.each do |vn, inf|
            if vn.sub(/\(.+\)\Z/,"") == v
              info = inf
              break
            end
          end
          unless info
            pp vars
            raise "parse error: #{group} #{v}"
          end
        end
        nm_params[v] ||= Array.new
        nm_params[v].push [modname, group]
        file.print <<EOL
!>        <tr>
!>          <td>#{v}</td>
!>          <td>#{info[:type]}</td>
!>          <td>#{info[:val]}</td>
!>          <td>#{info[:comments].join('<br />')}</td>
!>        </tr>
EOL
      end
      file.print "!>      </table><br>\n"
    end

    file.print <<EOL
!>    </ul>
!
!> @anchor history_#{modname}
!> @par History Output
EOL
    if hist.empty?
      file.print "!> No history output\n"
    else
      file.print <<EOL
!>    <table border=1>
!>      <tr><th>name</th><th>description</th><th>unit</th><th>variable</th></tr>
EOL

    hist.sort.each do |name, info|
      name = get_data(name,data,vars)
      desc = get_data(info[:desc],data,vars)
      if /\A(.*)trim\(([^()]+)(?:\([^()]+\))?\)(.*)\Z/ =~ name || /\A()(TRACER_NAME)\([^()]+\)()\Z/ =~ name
        pre = $1
        vn = $2
        post = $3
        pre.sub!(/'?([^']+)'\/\//, '\1')
        post.sub!(/\/\/'([^']+)'?/, '\1')
        vn = "@ref #{vn.downcase} \"#{vn}\"" if vn == "TRACER_NAME"
        name = "#{pre}<em>{#{vn}}</em>#{post}"
      end
      if /\A(.*)trim\(([^()]+)(?:\([^()]+\))?\)(.*)\Z/ =~ desc || /\A()(TRACER_DESC)\([^()]+\)()\Z/ =~ desc
        pre = $1
        vn = $2
        post = $3
        pre.sub!(/'?([^']+)'\/\//, '\1')
        post.sub!(/\/\/'([^']+)'?/,'\1')
        vn = "TRACER_NAME" if vn == "TRACER_DESC"
        vn = "HYD_NAME" if vn == "HYD_DESC"
        vn = "AE_NAME" if vn == "AE_DESC"
        case vn
        when "HYD_NAME"
          if /NUM_NAME/ =~ name
            desc = "#{pre}<em>{#{vn}}</em>#{post};<br/><em>{NUM_NAME}</em> is NC, NR, NI, NS, NG, NH."
          else
            desc = "#{pre}<em>{#{vn}}</em>#{post};<br/><em>{#{vn}}</em> is QC, QR, QI, QS, QG, QH."
          end
        when "AE_NAME"
          desc = "#{pre}<em>{#{vn}}</em>#{post};<br/><em>{#{vn}}</em> is AD, ASO, AVA, AS, AR, ASS, AU, AT, AOC."
        else
          vn = "@ref #{vn.downcase} \"#{vn}\"" if vn == "TRACER_NAME"
          desc = "#{pre}<em>{#{vn}}</em>#{post};<br/><em>{#{vn}}</em> depends on the physics schemes, e.g., QV, QC, QR."
        end
      end
      if /\A([^()]+)\([^()]+\)%(.+)\Z/ =~ name 
        vn  = $1
        elm = $2
        ary = data[vn].map do |d|
          if /\AVinfo\((.+)\)\Z/ =~ d
            body = $1.split(",").map{|s| s.strip.sub!(/\A'(.+)'\Z/,'\1')}
            [ body[0], modname, body[1], body[2], body[0]]
          else
            raise d
          end
        end
      else
        unit = get_data(info[:unit],data,vars)
        var = get_data(info[:var],data,vars)
        if /\A([^()]+)\(.+\)\Z/ =~ unit
          unit = "<em>{#{$1}}</em>"
        end
        if /\A([^()]+)\(.+\)\Z/ =~ var
          var = "<em>{#{$1}}</em>"
        end
        ary = [
          [name, modname, desc, unit, var]
        ]
      end
      ary.each do |ar|
        history[ar[0]] ||= [ar[2], Array.new]
        history[ar[0]][1].push ar[1]
        file.print <<EOL
!>      <tr>
!>        <td>#{ar[0]}</td>
!>        <td>#{ar[2]}</td>
!>        <td>#{ar[3]}</td>
!>        <td>#{ar[4]}</td>
!>      </tr>
EOL
      end
    end
    file.print <<EOL
!>    </table>
EOL
    end
    file.print <<EOL
module #{modname}
end module
EOL
  end
end


system("mkdir -p #{output_dir}")

if nm_params.any?
  File.open("#{output_dir}/namelist.dox","w") do |file|
    file.print <<EOL
!> @page namelist NAMELIST Parameters
!> @section #{title} #{title.upcase}
!> <p>
!> Click the corresponding "NAMELIST Group" or "module name" to see detail information of each variable
!> </p>
!> <table>
!> <tr><th>Variable name</th><th>NAMELIST Group</th><th>module name</th></tr>
EOL
    nm_params.sort.each do |name,ary|
      list = ary.map do |mod, group|
        file.print "!>    <tr><td>#{name}</td><td>@ref namelist_#{mod}_#{group} \"#{group}\"</td><td>@ref #{mod.downcase} \"#{mod}\"</td></tr>\n"
      end
    end
    file.print <<EOL
!>    </table>
EOL
  end
end


if history.any?
  File.open("#{output_dir}/history.dox","w") do |file|
    file.print <<EOL
!> @page history History Variables
!> @section #{title} #{title.upcase}
!> <p>
!> Click the corresponding "module name" to see detail information of each variable.
!> </p>
!> <p>
!> Note: "(w/ XXXXX)" means that it is enabled only when the SCALE is compiled with "-DXXXXX" compile option.
!> </p>
!> <table>
!> <tr><th>Variable name</th><th>Description</th><th>module name</th></tr>
EOL
    history.sort.each do |name,dat|
      desc,list = dat
      list = list.map{|mod|
        "@ref history_#{mod} \"#{mod}\""
      }
      file.print "!>    <tr><td>#{name}</td><td>#{desc}</td><td>#{list.join(", ")}</td></tr>\n"
    end
    file.print <<EOL
!>    </table>
EOL
    if title == "scale-rm"
    file.print <<EOL
!> <p>
!> @section tracer_name TRACER_NAME
!> <table>
!>  <tr><th>Category</th><th>Scheme name</th><th>Tracers</th><th>File</th></tr>
!>  <tr><td rowspan="3">Cloud microphysics</td><td>Kessler</td><td>QV,QC,QR</td><td>scale_atmos_phy_mp_kessler</td></tr>
!>  <tr><td>Tomita08</td><td>QV,QC,QR,QI,QS,QG</td><td>scale_atmos_phy_mp_tomita08</td></tr>
!>  <tr><td>SN14</td><td>QV,QC,QR,QI,QS,QG,NC,NR,NI,NS,NG</td><td>scale_atmos_phy_mp_sn14</td></tr>
!>  <tr><td rowspan="1">PBL</td><td>MYNN</td><td>TKE_MYNN</td><td>scale_atmos_phy_bl_mynn</td></tr>
!>  <tr><td rowspan="1">SGS trublence</td><td>D1980</td><td>TKE_D1980</td><td>scale_atmos_phy_tb_d1980</td></tr>
!> </table>
!> </p>
EOL
    end
  end
end
