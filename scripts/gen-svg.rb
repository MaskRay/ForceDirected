#!/usr/bin/env ruby
# ./gen.rb <<EOF
# 4 5
# 0 1
# 1 2
# 0 2
# 0 3
# 0 4
# EOF
puts '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="440px" height="440px">'

f = ->x { gets.split.map &x }
n, m = f[:to_i]
es = m.times.map { f[:to_i] }
coords = nil

IO.popen '../src/force', 'r+' do |io|
  io.puts "#{n} #{m}"
  es.each {|e| io.puts e.join ' ' }
  io.close_write
  coords = n.times.map { io.gets.split.map &:to_f }
end

n.times do |i|
  x, y = coords[i]
  puts %Q{<circle cx="#{x+20}" cy="#{y+20}" r="5" fill="black"/>}
  puts %Q{<text x="#{x+25}" y="#{y+20}" fill="red">#{i}</text>}
end

m.times do |i|
  x, y = coords[es[i][0]]
  xx, yy = coords[es[i][1]]
  puts %Q{<line x1="#{x+20}" y1="#{y+20}" x2="#{xx+20}" y2="#{yy+20}" stroke="black"/>}
end

puts '</svg>'
