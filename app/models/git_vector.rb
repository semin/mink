class GitVector < ActiveRecord::Base

  include Git::Distance

  belongs_to  :scop_domain,
              :class_name   => "ScopDomain",
              :foreign_key  => "scop_id"

  has_one   :norm_git_vector

  has_many  :mink_git_similarities

  acts_as_network :similar_git_vectors,
                  :through                  => :git_vector_similarities,
                  :foreign_key              => 'git_vector_id',
                  :association_foreign_key  => 'similar_git_vector_id'

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
