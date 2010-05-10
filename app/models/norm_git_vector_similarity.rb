class NormGitVectorSimilarity < ActiveRecord::Base

  belongs_to  :norm_git_vector,
              :class_name   => "NormGitVector",
              :foreign_key  => "norm_git_vector_id"

  belongs_to  :norm_git_vector_target,
              :class_name   => "NormGitVector",
              :foreign_key  => "similar_norm_git_vector_id"

end
