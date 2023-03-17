import configparser



if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")
    headlines = [headlines for headlines in config["worklist_headlines"]]
    print(headlines)



