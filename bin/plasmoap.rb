#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'optparse'

$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'bio-plasmoap'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <fasta_file>

    Predict whether protein(s) are targeted to the apicoplast
    using the PlasmoAP algorithm described here:

    Dissecting apicoplast targeting in the malaria parasite Plasmodium falciparum.
    Foth BJ, Ralph SA, Tonkin CJ, Struck NS, Fraunholz M, Roos DS, Cowman AF, McFadden GI.
    Science. 2003 Jan 31;299(5607):705-8.
    PMID: 12560551

    It uses SignalP version 3 (not 4!), which it expects to be on the PATH.
\n"
end; o.parse!

# print out a list of proteins with yes/no answers
puts [
  'Name',
  'PlasmoAP Score',
  'Apicoplast Targeted',
  'Points'
].join("\t")

runner = Bio::PlasmoAP.new

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
