from pick import pick

title = 'Please choose your favorite programming language (press SPACE to mark, ENTER to continue): '
options = ['Java', 'JavaScript', 'Python', 'PHP', 'C++', 'Erlang', 'Haskell']
selected = pick(options, title, multiselect=True, min_selection_count=1)
print(selected[0][0])