#!/usr/bin/env ruby

require "pp"
require "./options_HiKinetics.rb"
require 'open3'


# if ARGV.length != 0 then
#   exec("./HiKinetics.rb -h")
#   exit(-1)
# end

# Cross-platform way of finding an executable in the $PATH.
#
#   which('ruby') #=> /usr/bin/ruby
def which(cmd)
  exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
  ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
    exts.each { |ext|
      exe = File.join(path, "#{cmd}#{ext}")
      return exe if File.executable? exe
    }
  end
  return nil
end


if (!which('RNAHeliCes'))
  puts ""
  puts "--- Could not detect the RNAHeliCes program which is needed in this script ---"
  puts ""
  puts "Please make sure that you have installed it before running this script."
  puts ""
  exit 1
end

if (!which('treekin_hi'))
  puts ""
  puts "--- Could not detect the treekin_hi program which is needed in this script ---"
  puts ""
  puts "Please make sure that you have installed it before running this script."
  puts ""
  exit 1
end

options = HiKineticsOptparser.parse(ARGV)
#pp options
inputfilename = options.inputfilename
if (inputfilename.nil?)
  puts ""
  puts "--- Please specify input file and output file rootname ---"
  puts ""
  puts "Please type ./HiKinetics.rb -h for help in usage."
  puts ""
  #puts "Could not detect the input file, please use option -i to read input. See -h for more possible usages."
  exit 1
end

outputrootname = options.outputrootname
if (outputrootname.nil?)
  outputrootname = inputfilename.split(".")[0]
end

if (options.sn) 
  x_thresh = 0
else
  x_thresh = 1000000
end
#puts outputrootname
command = Thread.new do
  #system(`curl #{url} -o /tmp/#{filename}`) # long-long programm
  puts "generating hishapes ..."
  system(`./0_generate_hishapes.pl #{inputfilename} #{outputrootname} #{options.kbest} #{options.hishape_type} #{x_thresh}`)
end
command.join                 # main programm waiting for thread

command = Thread.new do
  puts "calculating partition function ..."
  system(`./1_add_partition_openchain_to_res.pl #{inputfilename} #{outputrootname} #{options.initial_hishape}`)
end
command.join                 # main programm waiting for thread

command = Thread.new do
  puts "preparing faa files for hipath calculation ..."
  system(`./2_generate_faa.rb #{outputrootname}`)
end
command.join                 # main programm waiting for thread

command = Thread.new do
  puts "calculating hipaths ..."
  system(`./3_generate_hipath.pl #{outputrootname}`)
end
command.join                 # main programm waiting for thread

command = Thread.new do
  puts "generate gnuplot file ..."
  system(`./4_hipath_generate_dat_plt.rb #{outputrootname}`)
end
command.join                 # main programm waiting for thread

# # command = Thread.new do
#   puts "generate #{outputrootname}.ps showing folding kinetics ..."
#   system(`gnuplot #{outputrootname}.plt`)
# # end
# # command.join                 # main programm waiting for thread

# puts "generate #{outputrootname}.pdf showing folding kinetics ..."
# exec("ps2pdf #{outputrootname}.ps")


# Specific options:
#     -i, --input FILE                 Read input from specified filename
#     -o, --output STRING              Specify root name for output
#     -k, --kbest [HISHAPE NUMBER]     Choosing the [HISHAPE NUMBER] best hishapes for folding kinetics
#     -t, --type [HISHAPE TYPE]        Specify abstraction type of the best hishapes for folding kinetics
#     -p, --population [HISHAPE]       Set initial population of state [HISHAPE] for folding kinetics
#     -s, --sn                         Specify a hishape abstraction type
