class GiVector < ActiveRecord::Base

  include Gi::Distance

  belongs_to  :scop_domain,
              :class_name   => "ScopDomain",
              :foreign_key  => "scop_id"

  has_one   :norm_gi_vector

  has_many  :mink_gi_similarities

  acts_as_network :similar_gi_vectors,
                  :through                  => :gi_vector_similarities,
                  :foreign_key              => 'gi_vector_id',
                  :association_foreign_key  => 'similar_gi_vector_id'

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
