library(ImageGP)
library(multcompView)

data <- sp_readTable("vegan.txt")

head(data)

rownames(data)

sp_diff_test(data, stat_value_variable="ACE",
                            stat_group_variable="Group")

a = sp_diff_test(data[1:12,], stat_value_variable="ACE",
                            stat_group_variable="Group")

a$data$y

sp_diff_test(data, stat_value_variable="chao1",
                            stat_group_variable="Group")


data = sp_diff_test(data, stat_value_variable="chao1",
                            stat_group_variable="Group")
data

split(data, data$Group)

c = lapply(split(data, data$Group), function (x) {max(x[["chao1"]])})

c[data$Group]

sp_diff_test_group_vector(data,
                                         stat_value_variable="ACE",
                                         stat_group_variable="Group",
           group_variable = "Site")

do.call(rbind,lapply(split(data, data$Group), function (x) {max(x[["chao1"]])}))

split(data, data$Site)

do.call(rbind, lapply(split(data, 1), sp_diff_test,
                      stat_value_variable="ACE",
                      stat_group_variable="Group"))
