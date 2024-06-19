load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
  name = "constants",
  hdrs = ["constants.h"],
)

cc_library(
  name = "algebra",
  hdrs = [
    "algebra_element.h",
    "group_element.h",
    "manifold_element.h",
  ],
  deps = ["@eigen"],
)

cc_library(
  name = "lie_group",
  hdrs = [
    "lie_group_element.h",
    "lie_algebra_element.h",
  ],
  deps = [":algebra"],
)

cc_test(
  name = "test_algebra_element",
  srcs = ["test_algebra_element.cc"],
  deps = [
    ":algebra",
    "@gtest//:gtest_main",
  ],
)

cc_test(
  name = "test_group_element",
  srcs = ["test_group_element.cc"],
  deps = [
    ":algebra",
    "@gtest//:gtest_main",
  ],
)

cc_test(
  name = "test_manifold_element",
  srcs = ["test_manifold_element.cc"],
  deps = [
    ":algebra",
    ":constants",
    "//utils:angles",
    "@gtest//:gtest_main",
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