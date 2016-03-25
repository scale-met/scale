=begin
Document generator for SCALE

This script generate list of namelist group and history.
The document is output in the HTML format.

Usage: ruby mkdoc.rb doc_dir
=end
require "pp"

def usage
  print "Usage: #$0 doc_dir\n"
  exit
end

output_dir = "./list"

docdir = ARGV.shift || usage
topdir = File.join(docdir, "..")
srcdir = File.join(topdir, "src")


def parse_line(line)
  line_org = line.dup
  line.sub!(/^\s*&/, "")
  if /^([^!]*)!(.*)$/ =~ line
    line = $1
    comment = $2
    comment.sub!(/^\s*</,"") if comment
    comment.strip! if comment
  else
    comment = nil
  end

  if /^(.+)&\s*$/ =~ line
    cont = true
    reg = $1
  else
    cont = false
    reg = line
    if /&/ =~ line && /["'][^"']*&[^"']*["']/ !~ line
      p line_org
      p line
      raise "pearse error"
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

files.flatten!.sort!

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
        when /^\s*(public|private)\s* ::/
          next
        when /^\s*module\s*([^\s]+)\s*$/
          modname = $1.strip
        when /^([^,]+).*(intent\([^\)]+\))?.*::([^=]+)(=.*)?$/
          next if $2
          type = $1
          name = $3
          val = $4
          type.strip!
          name = name.strip.upcase
          val.sub!(/\A=/,"").strip! if val
          if /,/ =~ name
            nas = name.split(",")
            nas.each do |na|
              na.strip!
              vars[na] = {:type => type, :val => nil, :comments => comments}
            end
            vars[nas[-1]][:val] = val if val
          else
            vars[name] = {:type => type, :val => val, :comments => comments}
          end
        when /NAMELIST\s*\/([^\)]+)\/(.+)$/
          group = $1
          lists = $2
          group.strip!
          lists = lists.split(",").map{|c| c.strip.upcase}
          namelists[group] = lists
        when /call HIST_in\s*\((.+)$/
          next if modname == "scale_history"
          str = $1.strip.sub(/\)\Z/,"").strip
          str.gsub!(/\(:[^)]*\)/,'')
          info = str.split(",").map{|c| c.strip.sub(/\A'(.*)'\Z/,'\1')}
          hist[info[1]] = {:unit => info[3], :desc => info[2], :var => info[0]}
        when /data\s+([^\s]+)\s+\/(.+)\/$/
          data[$1] = $2.split(',').map{|s| s.strip.sub(/^'/,"").sub(/'$/,"")}
        end

      end # while line
    end # open
  end

  vars.update(Index){|k,v1,v2| v1}

  next if namelists.empty? && hist.empty?

  parent = File.dirname(mod_f)
  system("mkdir -p #{output_dir}/#{parent}")

  tree[parent] ||= Array.new
  tree[parent].push File.basename(mod_f)

  File.open("#{output_dir}/#{mod_f}.html","w") do |file|

    file.print <<EOL
<html>
  <head>
    <title>#{modname} module</title>
  </head>
 <body>
  <h1>#{modname} module</h1>
  <h2>NAMELIST</h2>
    <ul>
EOL
    namelists.each do |group, lists|
      file.print <<EOL
      <li id="#{group}">#{group}</li>
      <table border=1>
        <tr><th>name</th><th>type</th><th>default value</th><th>comment</th></tr>
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
        nm_params[v].push [mod_f, group]
        file.print <<EOL
        <tr>
          <td>#{v}</td>
          <td>#{info[:type]}</td>
          <td>#{info[:val]}</td>
          <td>#{info[:comments].join('<br />')}</td>
        </tr>
EOL
      end
      file.print "      </table><br>\n"
    end

    file.print <<EOL
    </ul>

  <h2>History</h2>
    <table border=1>
      <tr><th>name</th><th>description</th><th>unit</th><th>variable</th></tr>
EOL

    hist.sort.each do |name, info|
      name = get_data(name,data,vars)
      history[name] ||= Array.new
      history[name].push mod_f
      file.print <<EOL
      <tr>
        <td id="#{name}">#{name}</td>
        <td>#{get_data(info[:desc],data,vars)}</td>
        <td>#{get_data(info[:unit],data,vars)}</td>
        <td>#{get_data(info[:var],data,vars)}</td>
      </tr>
EOL
    end
    file.print <<EOL
    </table>
  </body>
</html>
EOL
  end
end


system("mkdir -p #{output_dir}")
File.open("#{output_dir}/index.html","w") do |file|
  file.print <<EOL
<html>
  <head>
    <title>SCALE Document Index</title>
  </head>
  <body>
    <h1>SCALE Document Index</h1>

    <h2>Table</h3>
    <ul>
      <li><a href="#module">Module list</a></li>
      <li><a href="#namelist">Namelist parameters</a></li>
      <li><a href="#history">History variables</a></li>
    </ul>


    <h2 id="module">Module List</h2>

    <ul>
EOL
  tree.each.sort.each do |parent, mods|
    file.print "      <li>#{parent}</li>\n"
    file.print "      <ul>\n"
    mods.each do |mod|
      file.print "        <li id=\"#{mod}\"><a href=\"./#{parent}/#{mod}.html\">#{mod}</a></li>\n"
    end
    file.print "      </ul>\n"
  end
  file.print <<EOL
    </ul>

    <h2 id="namelist">NAMELIST Parameters</h2>
    <ul>
EOL
  nm_params.sort.each do |name,ary|
    list = ary.map do |mod, group|
      "<a href=\"#{mod}.html\##{group}\">#{group} in #{File.basename(mod)}</a>"
    end
    file.print "    <li>#{name}: #{list.join(", ")}</li>\n"
  end
  file.print <<EOL
    </ul>
    <h2 id="history">History Variables</h2>
    <ul>
EOL
  history.sort.each do |name,ary|
    list = ary.map do |mod|
      "<a href=\"#{mod}.html\##{name}\">#{File.basename(mod)}</a>"
    end
    file.print "    <li>#{name}: #{list.join(", ")}</li>\n"
  end
  file.print <<EOL
    </ul>

  </body>
</html>
EOL
end
