#include "gtest/gtest.h"
#include "lie/group_element.h"

namespace mana {
// The real numbers under addition form a group. Implement this group for
// testing.
class RealNumber : public GroupElement<RealNumber> {
 public:
  RealNumber(double value) : value_(value) {}

  // Access underlying value.
  double Value() const { return value_; }

  // Implement GroupElement interface.
  static RealNumber IdentityImpl() { return RealNumber(0); }
  RealNumber InverseImpl() const { return RealNumber(-value_); }
  RealNumber ComposeImpl(const RealNumber& rhs) const {
    return RealNumber(value_ + rhs.value_);
  }

  double value_;
};

TEST(GroupElement, Identity) {
  RealNumber a = RealNumber::Identity();
  EXPECT_EQ(a.Value(), 0);
}

TEST(GroupElement, Inverse) {
  RealNumber a(3);
  RealNumber b = a.Inverse();
  EXPECT_EQ(b.Value(), -3);
}

TEST(GroupElement, Compose) {
  RealNumber a(3);
  RealNumber b(5);
  RealNumber c = a.Compose(b);
  EXPECT_EQ(c.Value(), 8);

  RealNumber d = b.Compose(a);
  EXPECT_EQ(d.Value(), 8);
}

}  // namespace mana