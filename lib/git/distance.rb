module Git
  module Distance

    module ClassMethods

      require "narray"

      def euclidean_distance_between(foo, bar)
        vec_foo = NVector[
          foo.measure1 ,
          foo.measure2 ,
          foo.measure3 ,
          foo.measure4 ,
          foo.measure5 ,
          foo.measure6 ,
          foo.measure7 ,
          foo.measure8 ,
          foo.measure9 ,
          foo.measure10,
          foo.measure11,
          foo.measure12,
          foo.measure13,
          foo.measure14,
          foo.measure15,
          foo.measure16,
          foo.measure17,
          foo.measure18,
          foo.measure19,
          foo.measure20,
          foo.measure21,
          foo.measure22,
          foo.measure23,
          foo.measure24,
          foo.measure25,
          foo.measure26,
          foo.measure27,
          foo.measure28,
          foo.measure29,
          foo.measure30
        ]

        vec_bar = NVector[
          bar.measure1 ,
          bar.measure2 ,
          bar.measure3 ,
          bar.measure4 ,
          bar.measure5 ,
          bar.measure6 ,
          bar.measure7 ,
          bar.measure8 ,
          bar.measure9 ,
          bar.measure10,
          bar.measure11,
          bar.measure12,
          bar.measure13,
          bar.measure14,
          bar.measure15,
          bar.measure16,
          bar.measure17,
          bar.measure18,
          bar.measure19,
          bar.measure20,
          bar.measure21,
          bar.measure22,
          bar.measure23,
          bar.measure24,
          bar.measure25,
          bar.measure26,
          bar.measure27,
          bar.measure28,
          bar.measure29,
          bar.measure30
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
