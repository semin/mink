class Scop < ActiveRecord::Base

  cattr_reader :version

  @@version = "1.75"

  acts_as_nested_set

  def to_param
    self.sunid
  end

  def self.factory_create!(opts={})
    case opts[:stype]
    when "root" then ScopRoot.create!(opts)
    when "cl"   then ScopClass.create!(opts)
    when "cf"   then ScopFold.create!(opts)
    when "sf"   then ScopSuperFamily.create!(opts)
    when "fa"   then ScopFamily.create!(opts)
    when "dm"   then ScopProtein.create!(opts)
    when "sp"   then ScopSpecies.create!(opts)
    when "px"   then ScopDomain.create!(opts)
    else raise "Unknown SCOP hierarchy: #{opts[:stype]}"
    end
  end

  def hierarchy
    case stype
    when "root" then "Root"
    when "cl"   then "Class"
    when "cf"   then "Fold"
    when "sf"   then "Superfamily"
    when "fa"   then "Family"
    when "dm"   then "Protein"
    when "sp"   then "Species"
    when "px"   then "Domain"
    else "Unknown"
    end
  end

  def scop_domains
    leaves.select(&:rpall)
  end

  def html_sunid_link
    "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sunid=#{sunid}"
  end

  def html_sccs_link
    "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs=#{sccs}"
  end

end # class Scop


class ScopRoot < Scop
end


class ScopClass < Scop
end


class ScopFold < Scop
end


class ScopSuperFamily < Scop
end


class ScopFamily < Scop
end


class ScopProtein < Scop
end


class ScopSpecies < Scop
end


class ScopDomain < Scop

  has_one :mink_vector,
          :foreign_key  => :scop_id

  has_one :norm_mink_vector,
          :foreign_key  => :scop_id

  has_many  :gi_vectors,
            :foreign_key  => :scop_id

  has_many  :norm_gi_vectors,
            :foreign_key  => :scop_id

  has_many  :git_vectors,
            :foreign_key  => :scop_id

  has_many  :norm_git_vectors,
            :foreign_key  => :scop_id

  named_scope :rep95, :conditions => { :rep95 => true }

  def self.find_all_by_pdb_code(pdb_code)
    find(:all, :conditions => ["sid like ?", "_#{pdb_code.downcase}%"])
  end

  def pdb_code
    sid[1..4].upcase
  end

  def ranges_on_chains
    # "2hz1 A:2-124, B:1-50" => [A:2-124, B:1-50]
    description.gsub(/^\S{4}\s+/, '').gsub(/\s+/, '').split(',')
  end

  def include?(residue)
    result = false
    ranges_on_chains.each do |range|
      raise "Empty description!" if range =~ /^\s*$/
      case range
      when /^(\S):$/ # F:
        chain_code = $1
        if residue.chain.chain_code == chain_code
          result = true
        end
      when /^-$/ # -
        true
      when /^(-?\d+)-(-?\d+)$/ # 496-581
        res_from  = $1.to_i
        res_to    = $2.to_i
        if ((res_from..res_to).include?(residue[:residue_code]))
          result = true
        end
      when /^(\S):(-?\d+)-(-?\d+)$/ # A:104-157
        chain_code  = $1
        res_from    = $2.to_i
        res_to      = $3.to_i
        if ((residue.chain[:chain_code] == chain_code) &&
            (res_from..res_to).include?(residue[:residue_code]))
          result = true
        end
      else
        raise "#{self.description} should be added to Scop class!"
      end # case
    end # each
    result
  end

  def scop_class
    parent.parent.parent.parent.parent.parent
  end

  def scop_fold
    parent.parent.parent.parent.parent
  end

  def scop_superfamily
    parent.parent.parent.parent
  end

  def scop_family
    parent.parent.parent
  end

  def scop_protein
    parent.parent
  end

  def scop_species
    parent
  end

  def structure
    residues.first.chain.model.structure
  end

  def big_image
    "/figures/scop/#{sid}-solo-500.png"
  end

  def small_image
    "/figures/scop/#{sid}-solo-100.png"
  end

end
