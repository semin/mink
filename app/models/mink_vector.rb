class MinkVector < ActiveRecord::Base

  include Mink::Distance

  belongs_to  :scop_domain,
              :class_name   => "ScopDomain",
              :foreign_key  => "scop_id"

  has_one :norm_mink_vector

  define_index do
    indexes sid
    indexes sunid
    indexes sccs
    indexes scop_class_description
    indexes scop_fold_description
    indexes scop_superfamily_description
    indexes scop_family_description
    indexes scop_protein_description
    indexes scop_species_description
    indexes scop_domain_description
  end

end
