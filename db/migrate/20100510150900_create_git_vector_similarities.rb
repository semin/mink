class CreateGitVectorSimilarities < ActiveRecord::Migration
  def self.up
    create_table :git_vector_similarities, :force => true do |t|
      t.belongs_to  :git_vector
      t.belongs_to  :similar_git_vector
      t.float       :distance
    end

    add_index :git_vector_similarities, :distance
    add_index :git_vector_similarities, [:git_vector_id, :similar_git_vector_id], :name => "git1_git2"
    add_index :git_vector_similarities, [:similar_git_vector_id, :git_vector_id], :name => "git2_git1"
  end

  def self.down
    drop_table :git_vector_similarities
  end
end
