server/app.r: app/app.r app/rds/*
	cp app/app.r server/
	cp -r app/rds/* server/rds/
