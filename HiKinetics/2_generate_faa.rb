#!/usr/bin/env ruby

require "pp"
require "open3"


if ARGV.length != 1 then
  STDERR.puts "Usage: #{$0} outputrootname"
  exit(-1)
end
# dirfilename = "files_to_be_processed.txt"


outputrootname = ARGV[0]
##############################################################
####  read filename from dirfile  #### 
# if (File.exist?(dirfilename))
#   begin 
#     dirfile = File.new(dirfilename, "r")    
#   rescue
#     STDERR.print "Could not open file #{dirfilename}!\n"
#     exit 1
#   end
# else
#   STDERR.print "File #{dirfilename} does not exist!\n"
#   exit 1
# end



# dirdata = dirfile.readlines() 
# dirdata.each do |seqfilename|
#     seqfilename.strip!
#     rootname=seqfilename.split(".")[0]
    output_pre2="#{outputrootname}.pre2"

    
################ DEL hited #################    
    if (File.exist?(output_pre2))
      begin 
	inputfile = File.new(output_pre2, "r")    
      rescue
	STDERR.print "Could not open file #{output_pre2}!\n"
	exit 1
      end
    else
      STDERR.print "File #{output_pre2} does not exist!\n"
      exit 1
    end
    
    begin
      resfile = File.new("#{outputrootname}.res", "w")
    rescue
      STDERR.print "Could not open file #{outputrootname}.res!\n"
      exit 1
    end

    ss_array = []
    inputdata = inputfile.readlines()
    
    seq = ""
    i=-1
    inputdata.each do |line|
      line.strip!
      
      if (line.length > 10) then
	resfile.puts(line)
	line_data = line.split(%r{\s+})
	
	if line_data.length != 1 and line_data[0] != "length" then  # line_data[1].to_i < min
	  #min = line_data[1].to_f
	  ss_array << line_data[0]
	elsif line_data.length == 1 then
	  seq = line_data[0]
	end
	i=i+1
      end
    end
    resfile.close
    
    i=0
    ss_array.each do|ss1|
      j=0
      ss_array.each do|ss2| 
	    begin
	      faafile = File.new("#{outputrootname}.faas/#{i}_#{j}.faa", "w")
	    rescue
	      STDERR.print "Could not open file #{outputrootname}.faas/#{i}_#{j}.faa!\n"
	      exit 1
	    end
	    faafile.puts("> ss_number #{i} vs #{j}")
	    faafile.puts(seq)
	    faafile.puts(ss1)
	    faafile.puts(ss2)
	    faafile.close
	j=j+1
      end
      i=i+1
    end
