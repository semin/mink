class CreateGiValues < ActiveRecord::Migration
  def self.up
    create_table :gi_values, :force => true do |t|
      t.float :min_length
      t.float :min_int12
      t.float :min_inta12
      t.float :min_int12_34
      t.float :min_inta12_34
      t.float :min_int12_a34
      t.float :min_inta12_a34
      t.float :min_int13_24
      t.float :min_inta13_24
      t.float :min_int13_a24
      t.float :min_inta13_a24
      t.float :min_int14_23
      t.float :min_inta14_23
      t.float :min_int14_a23
      t.float :min_inta14_a23
      t.float :min_int12_34_56
      t.float :min_int12_35_46
      t.float :min_int12_36_45
      t.float :min_int13_24_56
      t.float :min_int13_25_46
      t.float :min_int13_26_45
      t.float :min_int14_23_56
      t.float :min_int14_25_36
      t.float :min_int14_26_35
      t.float :min_int15_23_46
      t.float :min_int15_24_36
      t.float :min_int15_26_34
      t.float :min_int16_23_45
      t.float :min_int16_24_35
      t.float :min_int16_25_34
      t.float :max_length
      t.float :max_int12
      t.float :max_inta12
      t.float :max_int12_34
      t.float :max_inta12_34
      t.float :max_int12_a34
      t.float :max_inta12_a34
      t.float :max_int13_24
      t.float :max_inta13_24
      t.float :max_int13_a24
      t.float :max_inta13_a24
      t.float :max_int14_23
      t.float :max_inta14_23
      t.float :max_int14_a23
      t.float :max_inta14_a23
      t.float :max_int12_34_56
      t.float :max_int12_35_46
      t.float :max_int12_36_45
      t.float :max_int13_24_56
      t.float :max_int13_25_46
      t.float :max_int13_26_45
      t.float :max_int14_23_56
      t.float :max_int14_25_36
      t.float :max_int14_26_35
      t.float :max_int15_23_46
      t.float :max_int15_24_36
      t.float :max_int15_26_34
      t.float :max_int16_23_45
      t.float :max_int16_24_35
      t.float :max_int16_25_34
      t.float :submax_length
      t.float :submax_int12
      t.float :submax_inta12
      t.float :submax_int12_34
      t.float :submax_inta12_34
      t.float :submax_int12_a34
      t.float :submax_inta12_a34
      t.float :submax_int13_24
      t.float :submax_inta13_24
      t.float :submax_int13_a24
      t.float :submax_inta13_a24
      t.float :submax_int14_23
      t.float :submax_inta14_23
      t.float :submax_int14_a23
      t.float :submax_inta14_a23
      t.float :submax_int12_34_56
      t.float :submax_int12_35_46
      t.float :submax_int12_36_45
      t.float :submax_int13_24_56
      t.float :submax_int13_25_46
      t.float :submax_int13_26_45
      t.float :submax_int14_23_56
      t.float :submax_int14_25_36
      t.float :submax_int14_26_35
      t.float :submax_int15_23_46
      t.float :submax_int15_24_36
      t.float :submax_int15_26_34
      t.float :submax_int16_23_45
      t.float :submax_int16_24_35
      t.float :submax_int16_25_34
    end

  end

  def self.down
    drop_table :gi_values
  end
end
