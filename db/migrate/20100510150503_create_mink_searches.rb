class CreateMinkSearches < ActiveRecord::Migration
  def self.up
    create_table :mink_searches, :force => true do |t|
      t.float       :cutoff
      t.string      :uuid
      t.timestamp   :started_at
      t.timestamp   :finished_at
      t.float       :elapsed_time
      t.string      :status
      t.float       :progress,     :default => 0
      t.float       :area_a
      t.float       :r_half_a
      t.float       :std_a
      t.float       :area_p
      t.float       :r_half_p
      t.float       :std_p
      t.float       :mean
      t.float       :std_mb
      t.float       :kurtosis
      t.float       :skewness
      t.float       :area_e
      t.float       :std_e
      t.float       :is
      t.float       :ori_area_a
      t.float       :ori_r_half_a
      t.float       :ori_std_a
      t.float       :ori_area_p
      t.float       :ori_r_half_p
      t.float       :ori_std_p
      t.float       :ori_mean
      t.float       :ori_std_mb
      t.float       :ori_kurtosis
      t.float       :ori_skewness
      t.float       :ori_area_e
      t.float       :ori_std_e
      t.float       :ori_is
      # columns for Paperclip plugin
      t.string      :pdb_file_name
      t.string      :pdb_content_type
      t.string      :pdb_file_size
      t.datetime    :pdb_updated_at
    end
  end

  def self.down
    drop_table mink_searches
  end
end
