class CreateNormGitVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :norm_git_vector_similarities, :force => true do |t|
      t.belongs_to  :norm_git_vector
      t.belongs_to  :similar_norm_git_vector
      t.float       :distance
    end

    add_index :norm_git_vector_similarities, :distance
    add_index :norm_git_vector_similarities, [:norm_git_vector_id, :similar_norm_git_vector_id], :name => "norm_git1_norm_git2"
    add_index :norm_git_vector_similarities, [:similar_norm_git_vector_id, :norm_git_vector_id], :name => "norm_git2_norm_git1"
  end

  def self.down
    drop_table :norm_git_vector_similarities
  end
end
