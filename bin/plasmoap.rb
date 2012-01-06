
require 'bio-plasmoap'

if $0 == __FILE__
  runner = Bio::PlasmoAP.new
  
  # print out a list of proteins with yes/no answers
  puts [
    'Name',
    'PlasmoAP Score',
    'Apicoplast Targeted',
    'Points'
  ].join("\t")
  
  Bio::FlatFile.auto(ARGF).each do |seq|
    result = runner.calculate_score(seq.seq)
    to_print = [seq.definition, result.to_s]
    if result.apicoplast_targeted?
      to_print.push 1
    else
      to_print.push 0
    end
    to_print.push result.points
    puts to_print.join("\t")
  end
end
