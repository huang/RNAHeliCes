require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'

class HiKineticsOptparser

  CODES = %w[iso-2022-jp shift_jis euc-jp utf8 binary]
  CODE_ALIASES = { "jis" => "iso-2022-jp", "sjis" => "shift_jis" }

  #
  # Return a structure describing the options.
  #
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    # We set default values here.
    options = OpenStruct.new
    #options.inputfiles = []
    options.inputfilename = nil
    options.outputrootname = nil
    options.kbest = 100
    options.hishape_type = 4
    options.initial_hishape = "[_]"
    options.sn = false
    #options.all_to_all = false
    #options.version = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: ./HiKinetics.rb -i infile -o outfile_rootname [options]\n"\
                    "Examples: ./HiKinetics.rb -i Input/xbix.seq -o Output/xbix -k 10 -t 2\n"\
                    "          ./HiKinetics.rb -i Input/xbix.seq -o Output/xbix_SN -k 10 -t 2 -s\n"\
                    "          ./HiKinetics.rb -i Input/xbix.seq -o Output/xbix_SN_P -k 10 -t 2 -s -p [7.5,15.5]\n"\
                    "          ./HiKinetics.rb -i Input/spliced_leader.seq -o Output/spliced_leader -k 100 -t 4\n"\
                    "          ./HiKinetics.rb -i Input/spliced_leader.seq -o Output/spliced_leader_SN -k 100 -t 4 -s\n"\
                    "          ./HiKinetics.rb -i Input/c_di_GMP_riboswitch.seq -o Output/c_di_GMP_riboswitch_SN_P -k 100 -t 4 -s -p [25.5]"
      #opts.banner = "Usage: HiKinetics.rb [options]\nExamples: ./HiKinetics.rb -i Input/spliced_leader.seq -o Output/spliced_leader -k 100 -t 2\n          ./HiKinetics.rb -i Input/spliced_leader.seq -o Output/spliced_leader_SN -k 100 -t 2 -s\n          ./HiKinetics.rb -i Input/c_di_GMP_riboswitch.seq -o Output/c_di_GMP_riboswitch_SN -k 100 -t 4 -s -p [25.5]"

      opts.separator ""
      opts.separator "Specific options:"

      
      # Mandatory argument.
#       opts.on("-i", "--input FILE", Array,
#               "Read input from specified file") do |input|
#         options.inputfiles << input
#       end
      opts.on("-i", "--input FILE",
              "Read input from specified filename") do |inputfilename|
        options.inputfilename = inputfilename
      end

      opts.on("-o", "--output STRING",
              "Specify root name for output") do |outputrootname|
        options.outputrootname = outputrootname
      end
      
      opts.on("-k", "--kbest [HISHAPE NUMBER]",
              "Choosing the [HISHAPE NUMBER] best hishapes for folding kinetics") do |kbest|
          options.kbest = kbest
      end
      
      opts.on("-t", "--type [HISHAPE TYPE]",
              "Specify abstraction type of the best hishapes for folding kinetics") do |hishape_type|
          options.hishape_type = hishape_type
      end 
      
      # Set initial population of state <int> to <double>
      opts.on("-p", "--population [HISHAPE]",
              "Set initial population of state [HISHAPE] for folding kinetics") do |initial_hishape|
          options.initial_hishape = initial_hishape
      end 

      # Boolean switch.      
      opts.on("-s", "--sn",
              "Restrict the folding space to strictly negative structure") do |sn|
          options.sn = sn
      end 
      
#       opts.on("-a", "--all", "Calculate folding pathways all-to-all") do |a|
#         options.all_to_all = a
#       end


      opts.separator ""
      opts.separator "Common options:"

      # No argument, shows at tail.  This will print an options summary.
      # Try it and see!
      opts.on_tail("-h", "--help", "Produce help message") do
        puts opts
        exit
      end

    end

    opts.parse!(args)
    options
  end  # parse()

end  # class OptparseExample

#options = OptparseExample.parse(ARGV)
#pp options
