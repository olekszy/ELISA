while True:
        try:
            row = z[z[1].str.contains(wave)].index.values.astype(int)[0]
        except IndexError:
            con = input("table not defined. Do you want to continue? y or n ")
            if con in ['y', 'Y', 'yes', 'Yes', 'YES']:
                exit()``