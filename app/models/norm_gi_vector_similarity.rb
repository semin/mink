class NormGiVectorSimilarity < ActiveRecord::Base

  belongs_to  :norm_gi_vector,
              :class_name   => "NormGiVector",
              :foreign_key  => "norm_gi_vector_id"

  belongs_to  :norm_gi_vector_target,
              :class_name   => "NormGiVector",
              :foreign_key  => "similar_norm_gi_vector_id"

end
