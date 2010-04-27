module Gi
  module Distance

    module ClassMethods

      require "narray"

      def euclidean_distance_between(foo, bar)
        vec_foo = NVector[
          foo.int12,
          foo.inta12,
          foo.int12_34,
          foo.inta12_34,
          foo.int12_a34,
          foo.inta12_a34,
          foo.int13_24,
          foo.inta13_24,
          foo.int13_a24,
          foo.inta13_a24,
          foo.int14_23,
          foo.inta14_23,
          foo.int14_a23,
          foo.inta14_a23,
          foo.int12_34_56,
          foo.int12_35_46,
          foo.int12_36_45,
          foo.int13_24_56,
          foo.int13_25_46,
          foo.int13_26_45,
          foo.int14_23_56,
          foo.int14_25_36,
          foo.int14_26_35,
          foo.int15_23_46,
          foo.int15_24_36,
          foo.int15_26_34,
          foo.int15_26_34,
          foo.int16_23_45,
          foo.int16_24_35,
          foo.int16_25_34
        ]

        vec_bar = NVector[
          bar.int12,
          bar.inta12,
          bar.int12_34,
          bar.inta12_34,
          bar.int12_a34,
          bar.inta12_a34,
          bar.int13_24,
          bar.inta13_24,
          bar.int13_a24,
          bar.inta13_a24,
          bar.int14_23,
          bar.inta14_23,
          bar.int14_a23,
          bar.inta14_a23,
          bar.int12_34_56,
          bar.int12_35_46,
          bar.int12_36_45,
          bar.int13_24_56,
          bar.int13_25_46,
          bar.int13_26_45,
          bar.int14_23_56,
          bar.int14_25_36,
          bar.int14_26_35,
          bar.int15_23_46,
          bar.int15_24_36,
          bar.int15_26_34,
          bar.int15_26_34,
          bar.int16_23_45,
          bar.int16_24_35,
          bar.int16_25_34
        ]

        NMath::sqrt((vec_foo - vec_bar)**2)
      end
    end

    def self.included(klass)
      klass.extend(ClassMethods)
    end

    def euclidean_distance_to(other)
      self.class.euclidean_distance_between(self, other)
    end
  end
end
