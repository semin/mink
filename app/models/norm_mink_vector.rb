class NormMinkVector < ActiveRecord::Base

  include Mink::Distance

  belongs_to  :mink_vector

  has_many :mink_search_hits

  has_many  :norm_mink_vector_similarities

  acts_as_network :similar_norm_mink_vectors,
                  :through                  => :norm_mink_vector_similarities,
                  :foreign_key              => "norm_mink_vector_id",
                  :association_foreign_key  => "similar_norm_mink_vector_id"

  def sorted_similar_norm_mink_vectors
    similar_norm_mink_vectors.sort_by { |other|
      euclidean_distance_to other
    }
  end

end

