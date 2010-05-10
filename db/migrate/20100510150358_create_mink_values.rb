class CreateMinkValues < ActiveRecord::Migration
  def self.up
    create_table :mink_values, :force => true do |t|
      t.float :min_area_a
      t.float :min_r_half_a
      t.float :min_std_a
      t.float :min_area_p
      t.float :min_r_half_p
      t.float :min_std_p
      t.float :min_mean
      t.float :min_std_mb
      t.float :min_kurtosis
      t.float :min_skewness
      t.float :min_area_e
      t.float :min_std_e
      t.float :min_is
      t.float :max_area_a
      t.float :max_r_half_a
      t.float :max_std_a
      t.float :max_area_p
      t.float :max_r_half_p
      t.float :max_std_p
      t.float :max_mean
      t.float :max_std_mb
      t.float :max_kurtosis
      t.float :max_skewness
      t.float :max_area_e
      t.float :max_std_e
      t.float :max_is
      t.float :submax_area_a
      t.float :submax_r_half_a
      t.float :submax_std_a
      t.float :submax_area_p
      t.float :submax_r_half_p
      t.float :submax_std_p
      t.float :submax_mean
      t.float :submax_std_mb
      t.float :submax_kurtosis
      t.float :submax_skewness
      t.float :submax_area_e
      t.float :submax_std_e
      t.float :submax_is
    end
  end

  def self.down
    drop_table :mink_values
  end
end
