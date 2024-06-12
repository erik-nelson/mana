load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
  name = "constants",
  hdrs = ["constants.h"],
)

cc_library(
  name = "lie_group",
  hdrs = [
    "lie_algebra_element.h",
    "lie_group_element.h",
    "traits.h",
  ],
  deps = [
    ":constants",
    "//base:common",
    "@eigen",
  ],
)

cc_test(
  name = "test_lie_group",
  srcs = ["test_lie_group.cc"],
  deps = [
    ":lie_group",
    "@gtest//:gtest_main",
  ],
)