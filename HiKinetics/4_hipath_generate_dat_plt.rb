#!/usr/bin/env ruby

require "pp"
require "open3"


if ARGV.length != 1 then
  STDERR.puts "Usage: #{$0} outputrootname"
  exit(-1)
end
# dirfilename = "files_to_be_processed.txt"


outputrootname = ARGV[0]


T = 37
kT = 0.00198717*(273.15 + T)
puts "kT=#{kT}"
##############################################################
####  read filename from dirfile  #### 
# dirdata = dirfile.readlines() 
# dirdata.each do |seqfilename|
#     seqfilename.strip!
#     rootname=seqfilename.split(".")[0]
# 
#     begin
#       hitedfile = File.new("#{rootname}.res.fas.hited", "r")
#     rescue
#       STDERR.print "Could not open file #{rootname}.res.fas.hited!\n"
#       exit 1
#     end

    hi_number=0
    #initstru_energies = []
    File.foreach("#{outputrootname}.res") {|line|
        line.strip!
        line_data = line.split(%r{\s+})
        if line_data.length>1 && line_data[1]!="=" then
            #initstru_energies << line_data[1].to_f
            hi_number=hi_number+1
        end
    }
    puts "hi_number=#{hi_number}"

    abs_hienergy = Array.new(hi_number) { Array.new(hi_number) {0.0} }
    ss_number = 1
    0.upto(hi_number-1) do |index_i|
      0.upto(hi_number-1) do |index_j|
	last_line = `tail -n 1 #{outputrootname}.faas/#{index_i}_#{index_j}.hip`
	last_line.strip!
        last_line_data = last_line.split(%r{\s+})
	_colon_line = `grep ' 0:' #{outputrootname}.faas/#{index_i}_#{index_j}.hip`
	_colon_line.strip!
        _colon_line_data = _colon_line.split(%r{\s+})
	#puts _colon_line
	#### NOTE alternative 1: if taking data from #{outputrootname}.res.line_data[1] (free energy of hishrep) calculated with turner2004, but the barrier energy (last_line_data[-2].to_f) calculated with turner1999 ==> inconsistent
	####      alternative 2: if taking data from #{outputrootname}.bar.line_data[2] (ensemble energy of hishape) calculated with turner2004, but the barrier energy (last_line_data[-2].to_f) calculated with turner1999 ==> inconsistent
	#abs_hienergy[index_i][index_j] = ((initstru_energies[index_i] + last_line_data[-2].to_f)*100).round / 100.0 # absoluteB = absoluteS + BarrierEnergy
	####      alternative 3: if taking data _colon_line_data[1].to_f (free energy of hishrep) calculated with turner1999, the barrier energy (last_line_data[-2].to_f) calculated with turner1999 too ==> consistent
        abs_hienergy[index_i][index_j] = _colon_line_data[1].to_f + last_line_data[-2].to_f  # absoluteB = absoluteS + BarrierEnergy
	ss_number = ss_number + 1
      end
    end
    #pp abs_hienergy


    ####### relax the table ########
    0.upto(hi_number-1) do |i|
      0.upto(i) do |j|
	if i==j then  # ENSURE values in the diagonal are 0.0
	    abs_hienergy[i][j]=0.0
	else
	    if (abs_hienergy[i][j]>abs_hienergy[j][i]) then
	      abs_hienergy[i][j]=abs_hienergy[j][i]
	    else  # abs_hienergy[j][i] is bigger
	      abs_hienergy[j][i]=abs_hienergy[i][j]
	    end
	end
      end
    end
    
   
    
    ##########################################################
    ########### HiTree ########### ./ShowTree -f rw1.tree.dat
    ##                                    ps2pdf rw1.tree.dat.ps
    ##########################################################
    ####### output .tree.dat file ######## 
    # >ggcc
    # GGGGGGCCCCCC
    #      0    ((((....))))   -7.20  [6.5]
    #      1    .((((...))))   -4.20  [7]
    #      2    ((((...)))).   -5.00  [6]
    # S         0         0       6.1       6.9
    # S         1       6.1         0       6.1
    # S         2       6.9       6.1         0
#     begin
#       treedatfile = File.new("#{rootname}.tree.dat", "w")
#     rescue
#       STDERR.print "Could not open file #{rootname}.tree.dat!\n"
#       exit 1
#     end

    i=-2
#     treedatfile.puts(">#{rootname}")
    hishapes = []
    File.foreach("#{outputrootname}.res") do |line|
        line.strip!
        line_data = line.split(%r{\s+})
        if line_data.length>1 && line_data[1]!="=" then
#             treedatfile.printf("%6d    ",i);
# 	    treedatfile.puts("#{line_data[0]}    #{line_data[1]}    #{line_data[2]}    #{((-kT*Math.log(line_data[3].to_f))*100).round / 100.0}")
	    hishapes[i] = line_data[2]
        elsif line_data.length==1 then
#             treedatfile.puts(line)
	end
	i=i+1
    end
#     0.upto(hi_number-1) do |i|
# 	    treedatfile.print("S");
# 	    treedatfile.printf("%10d",i);
# 	    0.upto(hi_number-1) do |j|
# 		treedatfile.printf("%10.4g",abs_hienergy[i][j]);
# 	    end
# 	    treedatfile.printf("\n");
#     end
    
    
    ##########################################################
    ########### Kinetic analysis based on Hishapes ###########
    ########### input: rw1.rat and rw1.bar
    ########### output: rw1.dat
    ########### using 'gnuplot rw1.plt' generate pdf which used the data from rw1.dat   
    ##########################################################
    ####### output .bar file ########
    # CUGCGGCUUUGGCUCUAGCC
    #      0    ....((((........))))    -4.52    [12.5]
    #      1    (((.(((....))).)))..    -3.95    [9.5]
    #      2    ..........(((....)))    -2.23    [15.5]
    #      3    ..((.......)).......    -0.55    [8]
    #      4    ....................    0.0    [_]
    #      5    ..((.............)).    0.1    [11]
    #      6    ....(((..((....)))))    0.36    [13.5]
    #      7    ..((((........)).)).    0.87    [10.5]
    #      8    ...((....)).........    2.07    [7.5]
    #      9    (((.((.......)))))..    1.86    [10]
    #     10    ....((.((.......))))    2.18    [13]
    #     11    ((........))........    2.0    [6.5]
    #     12    ...((....))((....)).    4.5    [7.5,15.5]
    #     13    .((......))((....)).    6.5    [6.5,15.5]
    hishape_fenergies = []
    id_fenergies = {}
    begin
      disfile = File.new("#{outputrootname}.bar", "w")
    rescue
      STDERR.print "Could not open file #{outputrootname}.bar!\n"
      exit 1
    end      
#    File.foreach("#{outputrootname}.seq") {|line| disfile.puts(line)}
    i=-2
#    disfile.puts(">#{outputrootname}")
    File.foreach("#{outputrootname}.res") do |line|
        line.strip!
        line_data = line.split(%r{\s+})
        if line_data.length>1 && line_data[1]!="=" then
            disfile.printf("%6d    ",i);
            #Z disfile.puts(line)
	    hishape_G = ((-kT*Math.log(line_data[3].to_f))*100).round / 100.0
	    #hishape_G = line_data[3].to_f
	    disfile.puts("#{line_data[0]}    #{hishape_G}    #{line_data[2]}")
	    #Z hishape_fenergies[i] = hishape_fenergy
	    id_fenergies[i] = hishape_G
        elsif line_data.length==1 then
            disfile.puts(line)
	end
	i=i+1
    end
    disfile.close
    #id_fenergies.sort_by {|k,v| v}.reverse
    sorted_id_fenergies = id_fenergies.sort{|a,b| a[1] <=> b[1]} 
    #pp sorted_id_fenergies
    
    
    ####### output.rat file ########
    begin
      ratesout = File.new("#{outputrootname}.rat", "w")
    rescue
      STDERR.print "Could not open file.rat!\n"
      exit 1
    end     
    0.upto(hi_number-1) do |i|
	    0.upto(hi_number-1) do |j|
		#ratesout.printf("%10.4g",abs_hienergy[i][j] - initstru_energies[i]);
	        ratesout.printf("%10.4g",abs_hienergy[i][j]);
	        #Z ratesout.printf("%10.4g", Math.exp(-(abs_hienergy[i][j]-hishape_fenergies[i])/kT));
	    end
	    ratesout.printf("\n");
    end
    ratesout.close
    
    
    ####### output .dat file ########
    system("treekin_hi --p0 #{hi_number}=1 --t0=0.001 --ratesfile=#{outputrootname}.rat -m H < #{outputrootname}.bar > #{outputrootname}.dat")
    #puts %x[#{cmd}]
# http://www.der-schnorz.de/2010/09/gnuplot-colors-presentations-papers-and-contrast/
linestyles = []    
linecolors=["red","blue","forest-green","magenta","gray","black","dark-red","orange"]    #,"royalblue","dark-orange"
1.upto(8) do |i|
  1.upto(8) do |j|
    linestyles << "set style line #{(i-1)*8+j} lt #{j} lc rgb \"#{linecolors[(j-1+i-1)%8]}\" lw 3"
  end
end
1.upto(8) do |i|
  1.upto(8) do |j|
    linestyles << "set style line #{(i-1+8)*8+j} lt #{j} lc rgb \"#{linecolors[(j-1+i-1)%8]}\" lw 1"
  end
end
#pp linestyles
# linestyles    
# "set style line 1 lt 1 lc rgb \"red\" lw 3",
# set style line 2 lt 2 lc rgb "blue" lw 3
# set style line 3 lt 3 lc rgb "forest-green" lw 3
# set style line 4 lt 4 lc rgb "magenta" lw 3
# set style line 5 lt 5 lc rgb "dark-orange" lw 3
# set style line 6 lt 6 lc rgb "royalblue" lw 3
# set style line 7 lt 7 lc rgb "black" lw 3
# set style line 8 lt 8 lc rgb "dark-red" lw 3
# set style line 9 lt 9 lc rgb "orange-red" lw 3
# set style line 10 lt 10 lc rgb "gray" lw 3
#pp hishapes


    ####### output .plt file ########    
    begin
      plt_file = File.new("#{outputrootname}.plt", "w")
    rescue
      STDERR.print "Could not open file #{outputrootname}.plt!\n"
      exit 1
    end  
    
    # if outputrootname is directory, split it and take the 2nd part and write it in plt_file
    outputrootname_data = outputrootname.split("/")
    outputrootname_last_part = outputrootname_data[-1] 
    
    #plt_file.puts("set title 'Hishape based kinetic analysis (#{outputrootname_last_part})'")
    #plt_file.puts("set terminal postscript eps enhanced color")  # dashed 16
    plt_file.puts("set terminal postscript enhanced eps color dashed 16")
    plt_file.puts("set xlabel \"{/Times=12 {/Symbol m}s}\"")
    plt_file.puts("set xrange [0.001:100000000]")
    plt_file.puts("set logscale x")
    plt_file.puts("set ylabel \"{/Times=12 Population density}\"")
    plt_file.puts("set yrange [0:1]")
    plt_file.puts("set term postscript enhanced font 'Time-roman,10'")
    plt_file.puts("set term postscript")
    #plt_file.puts("set output '#{outputrootname_last_part}.pdf'")
    plt_file.puts("set output '#{outputrootname_last_part}.ps'")
    #plt_file.puts("set key right top spacing 1.4") #title 'Legend' box 1
    #plt_file.puts("set key left bottom Left title 'Legend' box 3")
    plt_file.puts("set key right top Left")
    plt_file.puts("set key width 3")
    plt_file.puts("set key height 5")

    #plt_file.puts("set key spacing 1.4")
    #plt_file.puts("set key Left")
    #plt_file.puts("set key reverse")
    #plt_file.puts("set key width 20")
    #plt_file.puts("set key height 15")

    
    0.upto(hishapes.size-1) do |i|
      if (hishapes[i]=="[\_]") then
	hishapes[i]="[\\_]"
      end
    end
#     0.upto(hishapes.size-1) do |i|
#       puts hishapes[i]
#     end    
    
    0.upto(linestyles.size-1) do |i|
      plt_file.puts linestyles[i]
    end
    plt_file.puts("plot '#{outputrootname_last_part}.dat' using 1:2 title '#{hishapes[0]}' with lines ls 1, \\")
    2.upto(hi_number-1) do |i|  # 
	plt_file.puts("'#{outputrootname_last_part}.dat' using 1:#{i+1} title '#{hishapes[i-1]}' with lines ls #{i}, \\")
    end
#     plt_file.puts("plot '#{outputrootname_last_part}.dat' using 1:#{sorted_id_fenergies[0][0]} title '#{hishapes[sorted_id_fenergies[0][0]]}' with lines linewidth 3, \\")
#     2.upto(8) do |i|  # hi_number-1
# 	plt_file.puts("'#{outputrootname_last_part}.dat' using 1:#{sorted_id_fenergies[i+1][0]} title '#{hishapes[sorted_id_fenergies[i-1][0]]}' with lines linewidth 3, \\")
#     end
    plt_file.puts("'#{outputrootname_last_part}.dat' using 1:#{hi_number+1} title '#{hishapes[hi_number-1]}' with lines ls #{hi_number}")  # #{hishapes[hi_number-1]}
    #plt_file.puts("'#{outputrootname_last_part}.dat' using 1:#{hi_number+1} title '[\\_]' with lines ls #{hi_number}")  # #{hishapes[hi_number-1]}

    #plt_file.puts("pause -1 'Hit any key to continue'")
#end   
