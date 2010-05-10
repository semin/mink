class NormGitVector < ActiveRecord::Base

  include Git::Distance

  belongs_to  :git_vector

  belongs_to  :scop_domain,
              :class_name   => "ScopDomain",
              :foreign_key  => "scop_id"

  has_many  :git_search_hits

  has_many  :norm_git_vector_similarities

  acts_as_network :similar_norm_git_vectors,
                  :through                  => :norm_git_vector_similarities,
                  :foreign_key              => "norm_git_vector_id",
                  :association_foreign_key  => "similar_norm_git_vector_id"

  def sorted_similar_norm_git_vectors
    similar_norm_git_vectors.sort_by { |other|
      euclidean_distance_to other
    }
  end

end

