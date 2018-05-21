# DB設計

## usersテーブル

|Column|Type|Options|
|------|----|-------|
|name|string|null: false, index: true|
|email|string|null: false, unique: true|

### Association
- has_many :posts, dependent: :destroy
- has_many :likes, dependent: :destroy
- has_many :like_posts, dependent: :destroy
- has_many :comments, dependent: :destroy
- has_many :tags
- has_many :follows, dependent: :destroy
- has_many :followers, dependent: :destroy
- has_many :groups, through: groups_users, dependent: :destroy
- has_many :group_users, dependent: :destroy


## postsテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null: false|
|image_name|string||
|user_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- has_one    :group
- has_many   :likes, dependent: :destroy
- has_many   :comments
- has_many   :tags


## groupsテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- has_many :users, through: group_users
- has_many :group_messages
- has_many :group_users


## group_usersテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :group


## group_messagesテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null :false|
|image_name|string||
|group_id|integer|null :false, index: true, foreign_key: true|
|user_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :group
- belongs_to :user


## likesテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :post

## like_postsテーブル

|Column|Type|Options|
|------|----|-------|
|like_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :like
- belongs_to :post


## commentsテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null :false|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :post


## tagsテーブル

|Column|Type|Options|
|------|----|-------|
|type_id|integer|null :false, index: true, foreign_key: true|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :type
- belongs_to :user
- belongs_to :post


## typesテーブル

|Column|Type|Options|
|------|----|-------|
|type|string|null :false, index: true|

### Association
- belongs_to :tag


## relationshipsテーブル

|Column|Type|Options|
|------|----|-------|
|follower_id|integer|null :false, index: true, foreign_key: true, unique: true|
|followed_id|integer|null :false, index: true, foreign_key: true, unique: true|

### Association
- belongs_to :user
















