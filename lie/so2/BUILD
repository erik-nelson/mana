load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")
package(default_visibility = ["//visibility:public"])

cc_library(
  name = "so2",
  hdrs = [
    "so2_group_element.h",
  ],
  srcs = [
    "so2_group_element.cc",
  ],
  deps = [
    "@eigen",
    "//lie/base:lie_group",
    "//utils:angles",
  ],
)