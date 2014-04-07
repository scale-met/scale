
topdir = File.join(File.dirname(__FILE__), "..")
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
    elsif /\.F90\Z/ =~ d
      files.push d
    end
  end
  files
end

files = parse_dir(srcdir)

files.flatten!.sort!

tree = Hash.new
files.each do |fname|

  mod_f = fname.sub(/#{srcdir}\//,"").sub(/\.F90\Z/,"")

  vars = Hash.new
  namelists = Hash.new
  hist = Hash.new
  modname = nil
  File.open(fname) do |file|
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
      when /call HIST_(?:reg|in)\s*\((.+)$/
        next if modname == "scale_history"
        str = $1.strip.sub(/\)\Z/,"").strip
        str.sub!(/,[^,]+\Z/,"")
        info = str.split(",").map{|c| c.strip.sub(/\A'(.*)'\Z/,'\1')}
        hist[info[-3]] = {:unit => info[-1], :desc => info[-2], :var => info[0..-4].join(",")}
      end
    end
  end

  next if namelists.empty? && hist.empty?

  parent = File.dirname(mod_f)
  system("mkdir -p html/#{parent}")

  tree[parent] ||= Array.new
  tree[parent].push File.basename(mod_f)

  File.open("html/#{mod_f}.html","w") do |file|

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
      <li>#{group}</li>
      <table>
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
            require "pp"
            pp vars
            raise "parse error: #{group} #{v}"
          end
        end
        file.print <<EOL
        <tr>
          <td>#{v}</td>
          <td>#{info[:type]}</td>
          <td>#{info[:val]}</td>
          <td>#{info[:comments].join('<br />')}</td>
        </tr>
EOL
      end
      file.print "      </table>\n"
    end

    file.print <<EOL
    </ul>

  <h2>History</h2>
    <table>
      <tr><th>name</th><th>description</th><th>unit</th><th>variable</th></tr>
EOL
    hist.sort.each do |name, info|
      file.print <<EOL
      <tr>
        <td>#{name}</td>
        <td>#{info[:desc]}</td>
        <td>#{info[:unit]}</td>
        <td>#{info[:var]}</td>
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


File.open("html/index.html","w") do |file|
  file.print <<EOL
<html>
  <head>
    <title>SCALE Document Index</title>
  </head>
  <body>
    <h1>SCALE Document Index</h1>
    <ul>
EOL
  tree.each.sort.each do |parent, mods|
    file.print "      <li>#{parent}</li>\n"
    file.print "      <ul>\n"
    mods.each do |mod|
      file.print "        <li><a href=\"./#{parent}/#{mod}.html\">#{mod}</a></li>\n"
    end
    file.print "      </ul>\n"
  end
  file.print <<EOL
    </ul>
  </body>
</html>
EOL
end
