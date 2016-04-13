#!/usr/bin/ruby
# coding: utf-8

require 'bio'
require 'bio-signalp'

module Bio
  # PlasmoAP is a program for predicting apicoplast targetting sequences
  # in Plasmodium falciparum (the most fatal causative agent of malaria).
  # This algorithm was designed based on what was outlined in the 
  # PlasmoAP journal article 
  #
  # Dissecting apicoplast targeting in the malaria parasite Plasmodium falciparum.
  # Foth BJ, Ralph SA, Tonkin CJ, Struck NS, Fraunholz M, Roos DS, Cowman AF, McFadden GI.
  # Science. 2003 Jan 31;299(5607):705-8.
  # PMID: 12560551
  class PlasmoAP
    # Calculate the PlasmoAP score for a sequence (a string of amino acids)
    # sequence - the amino acids to test on
    # has_signal_sequence - Define if it has a signal sequence or not. The default
    # nil specifies that it should be worked out by SignalP.
    def calculate_score(sequence, has_signal_sequence = nil, signalp_cleaved_sequence = nil)
      # Only calculate signal sequence if it isn't already set by the parameter
      if has_signal_sequence.nil?
        # to_s means the sequence can be amino acid string or proper Bio::Sequence::AA object
        signalp = Bio::SignalP::Wrapper.new.calculate(sequence.to_s)
        unless signalp.kind_of?(Bio::SignalP::Version3::Result)
          raise "PlasmoAP uses signalp version 3, but different version (#{signalp.class}) was found"
        end
        
        has_signal_sequence = signalp.classical_signal_sequence?
        
        signalp_cleaved_sequence = signalp.cleave(sequence)
      elsif signalp_cleaved_sequence.nil?
        raise ArgumentError, "if the has_signal_sequence parameter is defined, then so must the signalp_cleaved_sequence be as well"
      end
      
      return PlasmoAPResult.new(0) if !has_signal_sequence #Both set of rules need a signal peptide
      
      cleaved = Bio::Sequence::AA.new(signalp_cleaved_sequence)
      
      set1 = set1?(cleaved)
      set2 = set2?(cleaved)
      additional = additional?(cleaved)
      
      points = 0
      points += 2 if set1
      points += 2 if set2
      points += 1 if additional
      return PlasmoAPResult.new(points)
    end
    
    private
    
    #      Set 1: a sequence is considered ‘positive’ if i) it starts with a signal peptide, ii) the 15
    #amino acids following the predicted signal peptide do not contain more than 2 acidic
    #residues, iii) the 80 amino acids following the predicted signal peptide cleavage site
    #contain a stretch of 40 amino acids with a total of at least 9 asparagines and/or lysines,
    #and iv) the asparagine/lysine-enriched region has a ratio of basic to acidic residues of at
    #least 5 to 3.
    def set1?(cleaved_amino_acids)
      set1 = false
      return false if cleaved_amino_acids.length <= 15
      aa = Bio::Sequence::AA.new(cleaved_amino_acids[0..14])
      if acidic_count(aa) <= 2 # ii)
        aa = Bio::Sequence::AA.new(cleaved_amino_acids[15..15+80-1])
        containing = nil
        # iii) contain a stretch of 40 amino acids with a total of at least 9 asparagines and/or lysines
        aa.window_search(40) do |window|
          if !containing #only the first one is needed
            comp = window.composition
            if comp['N']+comp['K'] >= 9
              containing = window
            end
          end
        end
        
        if containing
          if basic_count(containing).to_f / acidic_count(containing).to_f >= 5.0/3.0 # iv)
            set1 = true
          end
        end
      end
      
      return set1
    end
    
    
    #        Set 2: a sequence is considered ‘positive’ if it i) starts with a signal peptide, ii) if the 22
    #amino acids following the predicted signal peptide cleavage site exhibit a ratio of basic to
    #acidic residues of at least 10 to 7, iii) if the 80 amino acids following the predicted signal
    #peptide cleavage site contain a stretch of 40 amino acids with a total of at least 9
    #asparagines and/or lysines, and iv) if the asparagine/lysine-enriched region has a ratio of
    #basic to acidic residues of at least 10 to 9.
    def set2?(cleaved_amino_acids)
      set2 = false
      return false if cleaved_amino_acids.length <= 21
      aa = Bio::Sequence::AA.new(cleaved_amino_acids[0..21])
      if basic_count(aa).to_f / acidic_count(aa).to_f >= 10.0/7.0 #ii)
        
        # iii) if the 80 amino acids following the predicted signal
        #peptide cleavage site contain a stretch of 40 amino acids with a total of at least 9
        #asparagines and/or lysines
        aa = Bio::Sequence::AA.new(cleaved_amino_acids[22..22+80-1])
        containing = nil
        aa.window_search(40) do |window|
          if !containing #only the first one is needed
            comp = window.composition
            if comp['N']+comp['K'] >= 9
              containing = window
            end
          end
        end
        
        if containing
          if basic_count(containing).to_f / acidic_count(containing).to_f >= 10.0/9.0 # iv)
            set2 = true
          end
        end
      end
      return set2
    end
    
    
    # Additional point
    # basic => nil for none, true for basic, false for acidic
    def additional?(cleaved_amino_acids)
      cleaved_amino_acids.window_search(1) do |aa|
        if basic_count(aa) == 1
          return true
        elsif acidic_count(aa) == 1
          return false
        end
      end
      return nil
    end
    
    private
    ACIDIC_AMINO_ACIDS = %w(D E)
    BASIC_AMINO_ACIDS = %w(H K R)
    
    # Return the number of bases considered acidic in this amino acid sequence in neutral conditions.
    # Amino acids considered acidic are D and E. Input is a Bio::Sequence::AA object
    def acidic_count(amino_acid_sequence)
      count = 0
      comp = amino_acid_sequence.composition
      ACIDIC_AMINO_ACIDS.each do |acidic|
        if comp[acidic]
          count += comp[acidic]
        end
      end
      return count
    end

    # Return the number of bases considered basic in this amino acid sequence in neutral conditions.
    # Amino acids considered acidic are H, K and R. Input is a Bio::Sequence::AA object
    def basic_count(amino_acid_sequence)
      count = 0
      comp = amino_acid_sequence.composition
      BASIC_AMINO_ACIDS.each do |basic|
        if comp[basic]
          count += comp[basic]
        end
      end
      return count
    end


  end # End class PlasmoAP
  
  
  
  class PlasmoAPResult
    attr_reader :points
    def initialize(points)
      @points = points
      raise Exception, "Bad PlasmoAP Score points: #{points}" if points < 0 or points > 5
    end
    
    def to_s
      case @points
      when 0..2
        return '-'
      when 3
        return '0'
      when 4
        return '+'
      when 5
        return '++'
      end
    end  
    
    def ==(another)
      @points == another.points
    end
    
    # '+' or '++' scores were taken as apicoplast targeted in the paper
    # does this result pass that test?
    def apicoplast_targeted?
      @points >= 4
    end
    
    alias_method :predicted?, :apicoplast_targeted?
    alias_method :signal?, :apicoplast_targeted?
  end
end # End module Bio
