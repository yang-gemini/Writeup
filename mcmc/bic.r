n = 100
total = 0.0
for (i in 1:n) {
   total <- total + 2 / (2 + i - 1)
}
print(total)