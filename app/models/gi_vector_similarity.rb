class GiVectorSimilarity < ActiveRecord::Base

  belongs_to  :gi_vector,
              :class_name   => "NormGiVector",
              :foreign_key  => "gi_vector_id"

  belongs_to  :gi_vector_target,
              :class_name   => "NormGiVector",
              :foreign_key  => "similar_gi_vector_id"

end
