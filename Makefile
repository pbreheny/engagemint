server/app.r: app.r img/* rds/*
	cp app/app.r server/
	cp -r app/img/* server/img/*
	cp -r app/rds/* server/rds/*
