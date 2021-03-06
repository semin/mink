module Mink
  module Distance

    module ClassMethods

      require "narray"

      def euclidean_distance_between(foo, bar)
        vec_foo = NVector[
          foo.area_a,
          foo.r_half_a,
          foo.std_a,
          foo.area_p,
          foo.r_half_p,
          foo.std_p,
          foo.mean,
          foo.std_mb,
          foo.kurtosis,
          foo.skewness,
          foo.area_e,
          foo.std_e,
          foo.is
        ]

        vec_bar = NVector[
          bar.area_a,
          bar.r_half_a,
          bar.std_a,
          bar.area_p,
          bar.r_half_p,
          bar.std_p,
          bar.mean,
          bar.std_mb,
          bar.kurtosis,
          bar.skewness,
          bar.area_e,
          bar.std_e,
          bar.is
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
