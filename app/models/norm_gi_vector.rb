class NormGiVector < ActiveRecord::Base

  include Gi::Distance

  belongs_to  :gi_vector

  belongs_to  :scop_domain,
              :class_name   => "ScopDomain",
              :foreign_key  => "scop_id"

  has_many  :gi_search_hits

  has_many  :norm_gi_vector_similarities

  acts_as_network :similar_norm_gi_vectors,
                  :through                  => :norm_gi_vector_similarities,
                  :foreign_key              => "norm_gi_vector_id",
                  :association_foreign_key  => "similar_norm_gi_vector_id"

  def sorted_similar_norm_gi_vectors
    similar_norm_gi_vectors.sort_by { |other|
      euclidean_distance_to other
    }
  end

end

