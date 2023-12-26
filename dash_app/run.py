import os
from app import app

if __name__ == "__main__":
	port = int(os.environ.get("PORT", 3000))
	app.run_server(port=port, debug=True)